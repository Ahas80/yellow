# Audited methods for the longitudinal urinary *E. coli* analysis

## How to use this document

This is the defensible methods description supported by the current scripts and generated files. It deliberately separates:

1. the analysis that was actually executed on the mixed canonical assembly set;
2. the Longcycler-only sensitivity analysis;
3. exploratory analyses that should not be presented as confirmatory evidence; and
4. upstream laboratory or software details that are not recorded in the available files and therefore still need to be supplied.

The key methodological position is:

> The executed genomic pipeline selected one QC-passing canonical assembly per eligible participant episode, preferring Longcycler and using Flye only when the preferred Longcycler assembly was unavailable or failed the implemented QC criteria. This produced 556 selected episode-level genomes. A separate 532-row Longcycler-only analysis was then used to assess whether the longitudinal interpretation was sensitive to mixing assemblers.

Longcycler-only is therefore a **sensitivity analysis**, not a retrospectively relabelled primary pipeline. Restricting to one assembler reduces one possible source of technical variation, but it does not demonstrate that Longcycler is universally superior or remove every source of measurement error.

## 1. Study structure and analysis units

The project is a longitudinal observational analysis of repeated urinary *Escherichia coli* isolates. Several units occur in the pipeline and must not be treated as interchangeable:

- A **clinical record** is a row entering the clinical status map.
- A **clinical episode** is a participant–timepoint observation after episode consolidation and curation.
- A **participant** is a resident who may contribute several episodes.
- An **assembly record** is one candidate genome assembly. The same episode can have both Longcycler and Flye candidates.
- A **selected canonical row** is the single QC-passing assembly chosen to represent an episode in the executed genomic analysis.
- A **VF feature** is a binary virulence-factor-database gene feature, not an episode or participant.
- A **transition** is an adjacent pair of retained episodes within one participant after temporal ordering.
- A **pairwise comparison** is an unordered comparison between two selected genomes from the same participant; it is not restricted to adjacent episodes.

These distinctions explain why the numbers do not form one simple funnel. For example, 227 VF features are columns measured across genomes, whereas 394 transitions are repeated within-participant comparisons.

## 2. Clinical data preparation and UTI classification

Clinical exports from the available study batches were imported and harmonised in `00a_load_clean_clinical.R`. Participant identifiers, timepoints, event labels, culture information, urine-collection method and symptom fields were standardised before episode-level classification. Manual curation rules and duplicate handling were then applied through the shared pipeline helpers.

Episodes were classified in `00b_classify_episodes.R` using the current primary definition. A UTI required both:

1. urine-culture support at a lower-bound threshold of at least \(10^3\) CFU/mL; and
2. compatible symptoms under a catheter-aware rule.

For non-indwelling collection, the implemented symptom rule accepted a local urinary symptom or flank pain accompanied by a systemic symptom. For indwelling-catheter episodes, systemic symptoms were required. Episodes not satisfying both the culture and symptom components were labelled `Not_UTI`, with a subgroup retaining the reason for non-classification.

`Not_UTI` must not be described as “healthy”, “no infection” or “no bacteria”. In the current primary-included data, these rows represent bacteriuria that did not meet the operational symptomatic-UTI definition. Uncertain clinical evidence may also fall into an indeterminate `Not_UTI` subgroup and should remain visible in sensitivity analyses.

The current clinical status map contains 585 records from 167 participants before primary exclusions. Two records are excluded by the primary curation rules, leaving 583 episodes from 166 participants: 18 UTI and 565 Not_UTI. These counts use the primary `UTI_Status` field, not the older three-level `Infection_Status` field.

## 3. Assembly metadata, QC and canonical selection

`00_make_assembly_metadata.r` links prepared FASTAs to study overview and batch metadata and expands available Longcycler and Flye alternatives into assembly-level records. The documented workflow therefore begins with prepared clinical exports, metadata and FASTAs; it is not a raw-read-to-result workflow.

Assembly QC in `12a_wgs_qc.R` computes genome size, contig count, N50 and GC percentage from each available FASTA. An assembly passes the implemented screen when:

- the FASTA can be read without an error;
- total assembly length is between 4 and 6 Mb;
- the assembly contains no more than 200 contigs; and
- N50 is at least 20 kb.

Completeness and contamination constants present elsewhere in the configuration are not applied by this script and must not be claimed as part of this QC screen.

FASTA discovery identified 1,303 candidate FASTA files. This reconciles to 1,299 metadata rows rather than representing simple attrition: 1,295 rows have a matched existing FASTA and four are expected/no-FASTA audit rows, while eight discovered FASTAs are absent from the metadata table. Primary clinical/genomics curation leaves 1,291 candidate assembly records across 579 episode keys. Of these, 1,211 assembly records pass the implemented QC screen across 556 episode keys; the remaining 80 fail the assembly-size rule. Canonical selection then retains one QC-passing assembly per episode key, leaving 556 selected genomes from 162 participants: 532 Longcycler and 24 Flye fallback assemblies. Its clinical composition is 17 UTI and 539 Not_UTI episode rows.

## 4. Virulence-factor screening

`02_gene_presence_analysis.R` screens each selected canonical assembly with ABRicate against VFDB, using minimum identity and coverage thresholds of 80% each. Detected genes are converted to an episode-by-feature binary presence/absence matrix. The current canonical matrix contains 556 episode rows and 227 detected VFDB-derived features.

The 227 features must not all be called curated UPEC genes. Any apparent gene “gain” or “loss” between episodes means a difference in detected binary presence/absence. It does not by itself prove biological acquisition or deletion because lineage replacement, assembly differences and screening error can produce the same observation.

The installed ABRicate version, VFDB release and database update date are not frozen in the present outputs and must be added before manuscript submission.

## 5. MLST lineage context

`06_MLST.R` integrates provider and local sequence-type results. A provider SeqSphere call is preferred when `PercGoodTargets >= 95`; a labelled local result is used only as a fallback when the provider result is missing. Missing typing remains missing.

In the executed mixed canonical set, 533 of 556 rows have a usable preferred ST call, covering 83 ST labels. In the Longcycler-selected sensitivity set, 514 of 532 episode rows link to a preferred ST call, covering 80 ST labels; 509 of those calls are provider QC-passing calls and five are local fallbacks.

These are **MLST calls linked to Longcycler-selected episode rows**, not necessarily MLST results calculated from Longcycler assemblies. Provider provenance commonly contains Flye or combined Flye/Longcycler input. ST is used as broad lineage context: matching STs are compatible with relatedness but cannot, on their own, establish that two genomes are the same strain. The exact SeqSphere scheme and version, provider export manifest and allele-QC provenance remain required.

`22_vf_build_analysis_dataset.R` is the explicit episode-level integration step. It joins selected-canonical VF calls, the primary clinical label, preferred MLST calls and selected-assembly metadata by participant–timepoint; checks duplicate and unmatched joins; and writes the 556-row `vf_analysis_ready.csv` dataset. It does not perform association testing.

## 6. Pairwise genomic comparison and cache provenance

`11_compare_strains.R`, supported by `11_compare_strains_helpers.R`, compares every unordered pair of selected QC-passing canonical genomes belonging to the same participant. With the present 556 selected genomes this gives 963 within-participant pair comparisons.

For each pair the script records:

- assembly-to-assembly MUMmer `dnadiff` average identity and `TotalSNPs`;
- Mash distance;
- ST agreement;
- VF-profile Jaccard similarity.

The `TotalSNPs` values in `pairwise_metrics.csv` are **assembly-to-assembly `dnadiff` SNP differences**. They are not Parsnp core-genome SNP distances and are not wgMLST allele distances.

The repaired workflow permits only endpoints where `selected_canonical == TRUE`, `QC_PASS == TRUE` and the current FASTA exists. It no longer falls back silently to another assembler. Each FASTA is fingerprinted with SHA-256. A cache key includes the ordered pair identifier, normalised FASTA paths and both content hashes. Every reusable `dnadiff` report has a JSON sidecar containing the paths, hashes, sizes, modification times, command, software version, generation time and report hash. A cached result is reused only when the entire provenance record matches; otherwise the comparison is regenerated. The legacy key-only cache is preserved for audit but is never read by the repaired workflow.

The project uses an operational threshold of `TotalSNPs <= 25` as strong same-strain support. This is a predefined project rule, not a universal biological cutoff. A suitable citation or dataset-specific calibration is still required. The broader composite `Classification` field also incorporates ST and accessory-profile information and must not replace the strict SNP-supported interpretation.

## 7. Separate core-genome and pangenome branches

The pipeline also contains separate whole-cohort branches that must not be conflated with script 11:

- `12b_core_snp.R` stages the 556 selected canonical genomes, runs Parsnp with a random reference and duplicate filtering, calculates a separate SNP-distance matrix with `snp-dists`, and constructs a neighbour-joining tree. This is the core-genome branch and is not an input to the pairwise `dnadiff` values.
- `12c_panaroo.R` links or generates Prokka GFF annotations and runs Panaroo with strict cleaning and invalid-gene removal. This produces the pangenome/accessory-genome branch.

The exact Parsnp, `snp-dists`, Prokka, Panaroo, MUMmer and Mash versions must be recorded in the final reproducibility manifest. The Parsnp random-reference choice and its stale-output exit behaviour also need explicit reporting.

## 8. Longitudinal reconstruction

`24_vf_longitudinal_dynamics.R` orders retained episode rows within participant using collection date when available and labelled fallback ordering otherwise. It compares adjacent retained episode rows, calculates VF-profile Jaccard similarity, records detected presence/absence changes and attaches pairwise genomic context.

The executed mixed-canonical dataset produces 394 adjacent retained-episode transitions from 144 participants. These comprise 370 Not_UTI-to-Not_UTI, 10 Not_UTI-to-UTI, 13 UTI-to-Not_UTI and one UTI-to-UTI transition. Of the 394 mixed transitions, 138 meet the operational `dnadiff TotalSNPs <= 25` rule.

`scripts/rebuild_longcycler_sensitivity.R` rebuilds the Longcycler-only analysis after applying the assembler restriction. It does not filter the mixed transition table, because removing a Flye episode can make the Longcycler observations on either side newly adjacent. The restricted set contains 532 episode rows from 161 participants: 16 UTI and 516 Not_UTI. Of the 371 rebuilt transitions from 139 participants, 362 were already adjacent in the mixed timeline and nine become adjacent only after restriction. In total, 140 of 371 meet the operational `dnadiff TotalSNPs <= 25` rule; five of the nine Not_UTI-to-UTI transitions meet that rule.

“Adjacent retained episodes” is the correct description. These comparisons can span unequal time intervals and do not imply continuous bacterial observation between sampling dates.

## 9. Cross-sectional and statistical analyses

The main inferential script, `14_genotype_phenotype_model.R`, uses the 556-row mixed canonical VF-ready dataset with 17 UTI outcomes. Features with prevalence between 5% and 95% enter an initial Fisher screen. Features with nominal `p < 0.10`, or the top 50 when fewer meet that rule, are passed to a logistic mixed model containing the feature, timepoint, batch and a participant random intercept. A standard logistic model is used as a fallback when mixed-model fitting errors.

These analyses are exploratory rather than confirmatory because:

- only 17 UTI rows are available;
- feature selection occurs before the final model family;
- false-discovery correction is applied to the selected model subset rather than every screened feature;
- singular mixed-model fits are retained; and
- the fallback model does not account for repeated participant observations.

No feature reaches FDR < 0.05 in the current mixed-model output. Six fits are singular and 37 are flagged for sparsity or separation. These results are hypothesis-generating and must not be presented as estimates of causal risk.

## 10. Claims outside the current evidence boundary

The claim-linked workflow does not reconstruct antibiotic exposures, demographic or host-factor effects, wgMLST allele distances, defensible mutation locations or a true AMR analysis. It therefore does not support claims that antibiotics caused strain replacement, that a detected VF caused UTI, that colonisation was protective, or that host factors drove a genomic transition. These are evidence boundaries rather than allegations about intent.

## 11. Reproducibility boundaries and missing information

The available scripts support reproducibility from prepared clinical exports, assembly metadata and pre-existing FASTAs onward. They do not document or execute the complete upstream laboratory-to-assembly workflow. Before the methods are finalised, the following must be supplied from laboratory records, run manifests or the original analyst:

- Oxford Nanopore platform, flow cell and library kit;
- basecaller, version, model and basecalling settings;
- read-level QC, filtering and per-isolate coverage;
- Longcycler and Flye versions, command lines, polishing and any manual assembly selection;
- ABRicate version and VFDB release/date;
- SeqSphere/wgMLST/MLST scheme name, version and provider manifest;
- Prokka, Panaroo, Parsnp, `snp-dists`, MUMmer and Mash versions;
- a frozen R package/environment manifest; and
- a citation or internal calibration for the operational <=25 `dnadiff`-SNP threshold.

## 12. Short manuscript-style version

Clinical records from the available study batches were harmonised and consolidated by participant and sampling episode. Episodes were classified as UTI when urine culture supported bacterial growth at a lower-bound threshold of at least \(10^3\) CFU/mL and compatible symptoms were present under a catheter-aware symptom rule; other bacteriuric episodes were classified as Not_UTI with reason-specific subgroups. After primary exclusions, 583 episodes from 166 participants were retained, comprising 18 UTI and 565 Not_UTI episodes.

FASTA discovery identified 1,303 candidate files. Reconciliation produced 1,299 metadata rows, followed by 1,291 primary candidate assembly records and 1,211 QC-passing records. Candidate Longcycler and Flye assemblies were evaluated using assembly length, contig count, N50 and FASTA readability. QC-passing assemblies were 4–6 Mb, contained no more than 200 contigs and had N50 >=20 kb. One canonical assembly was selected per participant–timepoint, preferentially retaining Longcycler when both assemblers passed. The executed mixed canonical analysis comprised 556 genomes from 162 participants, including 532 Longcycler and 24 Flye fallback assemblies. A Longcycler-only restriction retaining 532 genomes from 161 participants was analysed separately as an assembler-sensitivity analysis.

Selected assemblies were screened with ABRicate/VFDB at >=80% identity and >=80% coverage to generate 227 binary VFDB-derived features. Provider SeqSphere MLST calls meeting >=95% good-target QC were preferred, with labelled local calls used only when provider typing was absent; ST was interpreted as lineage context rather than strain identity. All unordered within-participant genome pairs were compared using assembly-to-assembly MUMmer `dnadiff` SNP differences, Mash distance, ST agreement and accessory-profile similarity. Current FASTAs were SHA-256 fingerprinted, and `dnadiff` results were reused only when path-, hash-, command- and report-level provenance matched. `TotalSNPs <=25` was treated as an operational threshold for strong same-strain support, pending external citation or dataset-specific calibration.

Longitudinal analyses ordered retained episodes by collection date where available and compared adjacent retained episodes within participant. The mixed canonical analysis contained 394 transitions from 144 participants, of which 138 met the operational <=25 `dnadiff`-SNP rule. Rebuilding after the Longcycler-only restriction produced 371 transitions from 139 participants, of which 140 met the rule; five of nine Not_UTI-to-UTI transitions met it. VF presence/absence differences were interpreted descriptively and did not, by themselves, establish biological gene acquisition, strain replacement or causation of UTI.

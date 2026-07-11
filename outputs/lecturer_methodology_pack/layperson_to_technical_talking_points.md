# rUTI methodology: layperson-to-technical talking points

## How to use this guide

This guide follows the eight core lecturer slides and the six appendix slides. Begin with the plain-language explanation under each slide. Add the technical layer only when the audience is ready or asks for more detail.

The central framing is:

> The analysis that has actually been run uses one quality-controlled, selected genome assembly per episode. This is the 556-row mixed canonical dataset. I also use a 532-row Longcycler-only restriction as a sensitivity analysis to ask whether mixing assemblers materially affects the interpretation.

Do not describe the Longcycler-only set as the completed primary pipeline. It is a single-assembler sensitivity analysis unless the whole project is formally rerun and redesignated.

## A 30-second explanation

> I start with urine-sampling episodes and determine whether each episode meets the study definition of urinary tract infection. I then link eligible episodes to bacterial genome assemblies, apply quality checks, and choose one assembly per episode. The executed genomic dataset contains 556 selected assemblies. Most are Longcycler assemblies, while 24 are Flye fallbacks. I analyse the mixed dataset and separately repeat key summaries using 532 Longcycler-selected episodes. The second analysis reduces one source of technical variation, but it does not imply that Flye is unusable or that Longcycler is universally best.

## A simple analogy

Imagine that sequencing reads are pieces of a very large jigsaw puzzle. An assembler is the set of instructions used to put the pieces together.

- Using the same instructions for every puzzle removes one possible reason for the finished pictures to differ.
- It does not make every puzzle identical or eliminate differences in the pieces, read quality, sequencing or later processing.
- Longcycler-only is therefore a useful consistency check, not a guarantee that every remaining difference is biological.

Safe wording:

> The single-assembler analysis removes one avoidable technical difference. It does not remove every source of uncertainty.

## Essential glossary

| Term | Plain-language meaning | Technical meaning in this project |
|---|---|---|
| Clinical record | One row supplied in the clinical data. | A recorded sampling event before the two primary exclusions. |
| Episode | One urine-sampling occasion with clinical context. | The main episode-level unit linked by participant and timepoint. |
| Participant | One person in the study. | A resident who may contribute several episodes. |
| Isolate | The bacterium selected from an episode. | A urinary *Escherichia coli* culture taken forward for genomic work. |
| Sequencing reads | Many short pieces of DNA information. | The raw sequence observations used as input to assembly. |
| Assembly | A reconstructed genome. | A FASTA sequence built from sequencing reads. FASTA is a standard text format for nucleotide sequences. |
| Assembler | The software route that reconstructs a genome. | Longcycler or Flye in the available assembly metadata. |
| Assembly record | One candidate reconstruction. | A row in the assembly/QC table; one episode can have alternative assembler records. |
| QC | Quality control. | Checks applied before deciding whether an assembly is eligible. |
| Canonical assembly | The one selected assembly for an episode. | A QC-passing record marked `selected_canonical == TRUE`. |
| WGS | Whole-genome sequencing. | Sequencing intended to represent the bacterial genome. |
| VF | Virulence factor. | A gene feature in the VFDB-derived presence/absence matrix. |
| VFDB | Virulence Factor Database. | The reference database screened by ABRicate. |
| MLST | Multilocus sequence typing. | A lineage-typing approach based on alleles at a small set of housekeeping loci. |
| ST | Sequence type. | The MLST lineage label, such as ST131 or ST73. It is not proof of strain identity. |
| SNP | Single-nucleotide polymorphism. | A one-base difference detected in a defined sequence comparison. |
| `dnadiff` | A genome-comparison program from MUMmer. | The tool that produces the assembly-to-assembly SNP differences in `pairwise_metrics.csv`. |
| Mash | A rapid genome-similarity method. | An additional distance estimate used alongside `dnadiff`. |
| Parsnp | A core-genome alignment program. | A separate branch used for a multi-genome core alignment and tree; it is not the source of `pairwise_metrics.csv`. |
| Panaroo | A bacterial pangenome program. | Software that builds gene-family presence/absence outputs from Prokka annotations. |
| wgMLST | Whole-genome multilocus sequence typing. | An allele-distance method with different units and error structure from SNP comparison. |
| CFU/mL | Colony-forming units per millilitre. | The study uses culture support at or above 10³ CFU/mL where quantitative data are available. |
| N50 | A measure of assembly continuity. | The contig length at which half the assembly is contained in contigs of that length or longer. |
| Transition | A before-and-after comparison. | Two episodes adjacent among the retained episodes for the same participant in pipeline time order. |
| Sensitivity analysis | A robustness check. | Repeating selected analyses after restricting the 556 mixed assemblies to 532 Longcycler assemblies. |

## Slide 1 — One episode becomes a genome profile, then profiles are compared over time

### Plain-language explanation

> The project links what was happening clinically at a urine-sampling episode to the genome of the *E. coli* isolated at that episode. When a participant has repeated episodes, I can compare the bacteria over time.

The presentation answers three questions:

1. How were UTI and Not_UTI episodes defined?
2. How was one usable genome selected for each eligible episode?
3. What can repeated bacterial genomes tell us about lineage, gene carriage and strain continuity?

### Technical layer

Clinical classification and genomic analysis are separate steps. The clinical label is created without using the genome result. Assembly choice then determines which genome representation enters genomic analysis.

The executed analysis contains 556 selected canonical assemblies. The Longcycler-only set contains 532 of these rows and is used as a sensitivity analysis.

### Safe wording

> Clinical status is defined first. Genomic data are then used to describe the bacteria and their relatedness, not to manufacture the clinical label.

### Avoid saying

- “The genome proves that an episode was a UTI.”
- “The Longcycler-only analysis is the only valid dataset.”
- “Any genomic difference must be biological.”

### Likely lecturer question

**What is the overall unit of analysis?**

> It changes by question. Clinical classification is episode-level, genomic features are assembly-linked episode rows, and longitudinal analyses use within-participant episode pairs. I state the unit beside every denominator.

## Slide 2 — Different numbers count different things

### Plain-language explanation

> The numbers do not form one simple shrinking funnel. Some count records, some count people, some count candidate assemblies, some count selected genomes, some count gene features and some count comparisons over time.

### Denominator ladder

| Count | Unit | Meaning |
|---:|---|---|
| 585 | Clinical records | Episodes classified before the two primary exclusions. |
| 583 | Included clinical episodes | Primary clinical denominator after exclusions. |
| 166 | Participants | People represented in the 583 included episodes. |
| 1,303 | Discovered candidate FASTAs | Genome files found during discovery; files are not independent episodes. |
| 1,299 | Assembly-metadata rows | Reconciliation table: 1,295 rows point to matched, existing FASTAs and four are audit-only rows with no FASTA. Eight of the 1,303 discovered FASTAs are absent from this table. |
| 1,291 | Primary assembly candidates | Longcycler/Flye alternatives retained for primary genomic curation, not 1,291 independent episodes. |
| 1,211 | QC-passing assembly records | Candidate records satisfying the implemented assembly-QC rules. |
| 556 | Selected canonical rows | One selected mixed-assembler genome representation per retained episode. |
| 162 | Participants | People represented in the 556-row mixed genomic set. |
| 532 | Longcycler-selected rows | Single-assembler sensitivity denominator. |
| 161 | Participants | People represented in the 532-row sensitivity set. |
| 227 | VFDB-derived features | Binary presence/absence columns, not samples. |
| 394 | Mixed longitudinal transitions | Adjacent retained episode pairs from 144 participants. |
| 371 | Longcycler-only transitions | Rebuilt adjacent pairs from 139 participants. |

### Technical layer

A participant can contribute several episode rows. An episode can contribute more than one candidate assembly record, but only one record is selected for the canonical genomic analysis. A participant with `n` retained episodes contributes `n - 1` adjacent transitions.

The 1,303-to-1,299 step is reconciliation, not a simple exclusion: 1,295 metadata rows refer to matched, existing FASTAs, four rows are retained only for audit visibility despite having no FASTA, and eight discovered FASTAs do not appear in the metadata table. The primary-analysis filter then leaves 1,291 candidate assembly rows, of which 1,211 pass the implemented QC checks.

### Safe wording

> Before interpreting a number, I ask what one row represents and what filter produced it.

### Likely lecturer question

**Why do the numbers change?**

> They change for two reasons: eligibility filters remove or retain rows, and the unit itself changes. For example, 227 is a number of genomic features, whereas 532 is a number of genome-linked episodes.

## Slide 3 — UTI requires culture support plus compatible symptoms

### Plain-language explanation

> An episode is labelled UTI only when the urine culture supports bacterial growth and the participant has compatible symptoms under the study’s catheter-aware rule. An episode that does not meet both requirements is labelled Not_UTI.

The clinical ladder is:

- 585 records were classified before primary exclusions.
- Two records were excluded: one secondary duplicate and one record with an unknown participant.
- 583 episodes from 166 participants remained.
- These comprised 18 UTI and 565 Not_UTI episodes.

### Important clarification

Not_UTI does not mean healthy, symptom-free, culture-negative or bacteria-free. In the current included dataset, the Not_UTI rows are bacteriuria episodes that do not meet the project’s complete UTI definition.

### Technical layer

Culture support uses a threshold of at least 10³ CFU/mL where quantitative CFU data are available. Compatible symptoms are interpreted using different logic for catheter and non-catheter episodes. Clinical status confidence and the component rules should be reported with the final methods.

The displayed 18/565, 17/539 and 16/516 splits use the current `UTI_Status` field. A legacy `Infection_Status` field is not the source of these presentation counts.

### Safe wording

> Not_UTI is the comparison category for bacteriuria episodes that did not satisfy the study’s combined culture-and-symptom UTI rule.

### Avoid saying

- “Not_UTI means no infection or no bacteria.”
- “The comparison group is healthy controls.”
- “The UTI definition was inferred from the genome.”

### Likely lecturer questions

**Why is the UTI count small?**

> The definition deliberately requires both culture support and compatible symptoms. This creates a narrow UTI group and a much larger, heterogeneous Not_UTI group.

**Why do I later see 17 or 16 UTI rows?**

> Eighteen is the included clinical denominator. Seventeen UTI rows have a selected canonical genome in the mixed analysis, and 16 remain after the Longcycler-only sensitivity restriction. The clinical definition has not changed; the genomic eligibility denominator has.

## Slide 4 — Clinical and assembly tracks converge before one genome is selected per episode

### Plain-language explanation

> The clinical track tells me which sampling episodes are eligible and whether each is UTI or Not_UTI. A separate assembly track finds and checks possible genome reconstructions. The two tracks are linked by participant and timepoint before one genome is selected for each retained episode.

### How to walk through the flowchart

Follow the numbered assembly/QC track. Stages 1–3 run left to right across the top row; the long diagonal arrow then carries the assembly track to stage 4, after which it continues left to right through stages 5 and 6:

1. **1,303 candidate FASTAs discovered.** These are files found on disk, not 1,303 people or episodes.
2. **1,299 assembly-metadata rows reconciled.** Of these, 1,295 point to matched, existing FASTAs and four are audit-only rows without a FASTA. Eight discovered FASTAs are absent from the metadata table, so this step is reconciliation rather than a simple four-file loss.
3. **1,291 primary assembly candidates retained.** Eight metadata rows do not meet the primary/genomics curation boundary.
4. **1,211 candidates pass implemented QC.** Eighty candidate records fail; in the current output these failures are due to the implemented genome-size rule.
5. **556 canonical genomes selected.** One eligible, QC-passing candidate is selected per retained episode; 655 other QC-passing alternatives are not selected because they represent additional candidate reconstructions rather than additional episodes.

Read the lower clinical lane in parallel:

- 583 included clinical episodes from 166 participants comprise 18 UTI and 565 Not_UTI episodes.
- Linking clinical episodes to selected canonical genomes yields 556 genomic episode rows from 162 participants: 17 UTI and 539 Not_UTI.
- The 27 included clinical episodes without a selected genome consist of four with no primary candidate and 23 with no QC-passing candidate.

Finally, the selected set splits into two reporting routes: the executed mixed analysis retains all 556 selected genomes, while the sensitivity route retains 532 selected Longcycler genomes, removes 24 Flye fallbacks and rebuilds the participant timelines.

### Implemented QC rules

An assembly passes the current implemented QC step when:

- the FASTA can be read;
- estimated genome size is between 4 and 6 megabases;
- it contains no more than 200 contigs; and
- N50 is at least 20 kilobases.

The primary assembly table contains 1,291 candidate records, of which 1,211 pass these implemented checks. Canonical selection produces 556 selected rows: 532 Longcycler and 24 Flye fallback rows.

### Technical layer

The headline steps shown on this slide are produced by `00a_load_clean_clinical.R` and `00b_classify_episodes.R` for the clinical track, and by `00_make_assembly_metadata.r` and `12a_wgs_qc.R` for assembly reconciliation, implemented QC and canonical selection. The project runner is not a raw-read-to-results workflow: it assumes that assembly FASTAs and assembly metadata already exist. Sequencing/read QC, Longcycler/Flye commands and polishing therefore remain upstream provenance that must be described separately.

The code contains completeness and contamination constants, but the current assembly-QC implementation does not apply those thresholds. Do not claim that these criteria were used unless they are added and rerun.

### Safe wording

> The downstream pipeline begins from available assembly FASTAs and metadata, applies the implemented contiguity and size checks, and selects one QC-passing canonical assembly per episode.

### Avoid saying

- “The R pipeline performed sequencing and assembly.”
- “All conventional genome-QC criteria were applied.”
- “There are 1,291 samples.”

### Likely lecturer questions

**Why is Flye used at all?**

> The canonical rule is Longcycler-preferred, with Flye used as a fallback where a selected Longcycler assembly is not available. The sensitivity analysis removes those 24 fallback rows to test whether the mixed route affects interpretation.

**Is Longcycler the best assembler?**

> The current analysis does not establish a universally best assembler. Longcycler-only is useful here because it standardises one component of the comparison. Assembler choice should also consider read type, assembly accuracy, structural resolution, polishing, downstream task and empirical validation.

## Slide 5 (REMOVED) — The executed analysis is mixed; Longcycler-only tests assembler sensitivity 

### Plain-language explanation

> The analysis that has already been run uses 556 selected assemblies. Most are Longcycler assemblies, with 24 Flye fallbacks. I compare that executed dataset with a 532-row Longcycler-only version to see whether the conclusions depend on mixing assembly methods.

| Executed mixed canonical analysis | Longcycler-only sensitivity analysis |
|---|---|
| 556 episode rows | 532 episode rows |
| 162 participants | 161 participants |
| 17 UTI and 539 Not_UTI rows | 16 UTI and 516 Not_UTI rows |
| 532 Longcycler plus 24 Flye fallbacks | Selected QC-passing Longcycler only |

The sensitivity analysis retains approximately 95.7% of the selected mixed rows.

### Technical layer

Restricting to one assembler can reduce assembler-associated heterogeneity, especially for small-distance genomic comparisons. It can also introduce selection if the excluded fallback rows differ systematically from retained rows. The sensitivity result must therefore be shown alongside the executed analysis, not presented as proof that exclusion is unbiased.

### Safe wording

> The mixed dataset is the executed analysis. Longcycler-only is a high-retention sensitivity analysis that removes one source of technical heterogeneity.

### Avoid saying

- “Longcycler-only is the clean truth.”
- “Flye rows are wrong.”
- “The sensitivity restriction removes all technical bias.”

### Likely lecturer question

**Am I losing too much information?**

> Twenty-four selected rows and one represented participant are removed. Most data are retained, but I will report both analyses because the excluded episodes may not be a random subset.

## Slide 6 — The selected genomes feed several parallel analysis branches

### Plain-language explanation

> The genome does not go through one single chain in which each result depends on the previous result. The selected assemblies feed several separate analyses, each answering a different question.

### Branch A: VFDB-derived feature screening

ABRicate screens selected assembly FASTAs against VFDB using 80% identity and 80% coverage thresholds. The analysis-ready matrix contains 227 binary detected/not-detected features.

Safe interpretation:

> These features describe detected VFDB-matched gene carriage. Presence does not establish expression, activity or causation.

### Branch B: MLST lineage context

Preferred MLST calls are linked to selected episode keys. In the Longcycler sensitivity set, 514 rows have a usable ST across 80 distinct ST labels. Of these typed rows, 509 use provider calls meeting the provider QC rule and five use labelled local fallback calls; 18 of the 532 sensitivity rows remain untyped.

The provider-call provenance often references Flye and/or Longcycler inputs. Therefore, the correct phrase is “Longcycler-selected episodes linked to preferred MLST calls,” not “MLST generated only from Longcycler assemblies.”

Safe interpretation:

> ST describes broad lineage context. The same ST is compatible with relatedness but does not prove that two isolates are the same strain.

### Branch C: pairwise assembly comparison

MUMmer `dnadiff` and Mash compare selected assemblies within the same participant. `pairwise_metrics.csv` contains assembly-to-assembly comparison metrics. These are not Parsnp core-genome SNP distances.

The project applies a predefined operational threshold of at most 25 `dnadiff` SNP differences as strong same-strain support. The threshold still needs an explicit citation or dataset-specific calibration and should not be described as universal.

### Integration bridge: one analysis-ready episode table

`22_vf_build_analysis_dataset.R` links the clinical status, selected assembly, VFDB-derived feature matrix and preferred MLST call at the episode level. This is an integration step, not an additional biological test: its purpose is to ensure that downstream summaries refer to the same participant/timepoint keys and denominators.

### Branch D: separate core-genome analysis

Parsnp aligns the multi-genome core using a selected reference, `snp-dists` creates distances from that core alignment, and a neighbour-joining tree provides population context. This branch must remain separate from `dnadiff` terminology and counts.

### Branch E: pangenome context

Prokka produces genome annotations and Panaroo constructs a strict-cleaned pangenome. This is a separate population/genome-content context branch. It is not the source of the pairwise `dnadiff` counts and does not by itself establish a causal mechanism.

The headline branch scripts are `02_gene_presence_analysis.R` for VFDB screening, `06_MLST.R` for lineage calls, `11_compare_strains.R` plus its helper for provenance-verified pairwise comparison, and `22_vf_build_analysis_dataset.R` for integration. `12b_core_snp.R` and `12c_panaroo.R` are named only where the separate core-genome and pangenome context is being explained.

### Likely lecturer questions

**Why is MLST not enough for same-strain inference?**

> MLST uses a small number of housekeeping loci. Many distinct strains can share a common ST. It is useful for lineage context, but strain-level inference needs higher-resolution evidence.

**Are SNP and wgMLST thresholds interchangeable?**

> No. SNP comparison counts nucleotide differences within a defined alignment/comparison, whereas wgMLST counts allele differences across a specific locus scheme. Their units, missing-data behaviour and technical error are different.

## Slide 7 (REMOVE)— Longitudinal comparisons must be rebuilt after changing retained episodes

### Plain-language explanation

> A transition compares one retained episode with the next retained episode from the same participant. If I remove a Flye episode from the middle of a timeline, two Longcycler episodes that were not previously neighbours can become adjacent. I therefore rebuild the timeline instead of merely filtering the old transition table.

### Verified transition counts

- The executed mixed analysis has 394 adjacent transitions from 144 participants.
- The Longcycler-only reconstruction has 371 transitions from 139 participants.
- Of these, 362 were already Longcycler-to-Longcycler adjacencies in the mixed timeline and nine are newly adjacent after Flye removal.
- The Longcycler reconstruction contains nine Not_UTI-to-UTI transitions.

### Current provenance-validated SNP results

- Mixed adjacent transitions meeting the operational threshold: 138 of 394.
- Longcycler-only adjacent transitions meeting the operational threshold: 140 of 371.
- Longcycler-only Not_UTI-to-UTI transitions meeting the operational threshold: five of nine.

These values were recomputed from the current selected assemblies using provenance-validated pairwise reports. The pairwise output contains 963 unordered within-participant comparisons from the 556 selected canonical genomes; every report is tied to the current endpoint FASTA paths, SHA-256 hashes and a matching provenance sidecar.

### Technical layer

“Adjacent” means adjacent among retained episodes in pipeline time order. It does not guarantee uninterrupted biological sampling, equal elapsed time or absence of unobserved intermediate strains.

Detected VF gains and losses are changes in presence/absence calls between retained episode rows. They are not direct proof that a gene was biologically acquired or deleted.

### Safe wording

> The result is consistent with genomic continuity between retained timepoints under the project’s operational threshold.

### Avoid saying

- “The same strain definitely persisted continuously between visits.”
- “Colonisation caused the later UTI.”
- “A detected gene gain proves horizontal acquisition.”

### Likely lecturer question

**What exactly does at most 25 SNPs mean?**

> It is the project’s operational interpretation boundary for the current `dnadiff` comparison. It provides strong support under this analysis framework, but it does not by itself prove transmission, direction, continuous carriage or causation.

**How can the Longcycler analysis have fewer total transitions but more threshold-supported transitions?**

> The Longcycler transition table is rebuilt, not obtained by taking a simple subset of the mixed table. Removing a Flye episode can make two surrounding Longcycler episodes newly adjacent, so the set of pairs being counted changes. The 140 and 138 values therefore come from two different adjacency structures and should be compared as sensitivity results, not treated as nested numerators.

## Slide 8 — The methods support description and qualified relatedness, not causation

### What is supported

- The clinical classification rule and its denominators.
- The implemented assembly-QC and canonical-selection process.
- VFDB-derived feature presence/absence summaries.
- Preferred MLST lineage context.
- Provenance-validated assembly-to-assembly relatedness summaries.
- Rebuilt within-participant adjacent transitions.
- A Longcycler-only sensitivity analysis alongside the executed mixed analysis.

### What is not supported by this workflow alone

- A claim that antibiotics caused strain replacement or UTI.
- A claim that colonisation was protective or caused later disease.
- A recalculated wgMLST result without the allele-distance table and scheme version.
- A demographic or host-factor conclusion without the required source data and repeated-measures analysis.
- An antibiotic-resistance conclusion without an appropriate resistance dataset and validated analysis.
- A claim that a particular sequence difference has a confirmed gene location without sequence-aware annotation and validation.

### Decisions to request from the lecturer

1. Confirm that the thesis should present the 556-row mixed canonical analysis as executed and Longcycler-only as sensitivity.
2. Agree how the operational 25-SNP boundary will be justified or calibrated.
3. Confirm which exploratory analyses belong in the thesis and which should remain supplementary.
4. Complete the missing upstream sequencing, assembly, database and software-version provenance.

### Closing wording

> The strongest defence is not that every result is definitive. It is that each result is tied to a clear input, filter, analysis unit and method, and that the interpretation stops where the available evidence stops.

## Appendix talking points

### Appendix A — Full denominator and filter table

For every number, show four items:

1. source file;
2. row filter;
3. unit and denominator; and
4. calculation or formula.

This prevents stale summaries from becoming sources of truth.

### Appendix B — QC and software thresholds

Distinguish:

- thresholds that the current code actually applies;
- constants that exist but are not applied;
- upstream procedures not represented in the R runner; and
- interpretation thresholds that still require scientific justification.

### Appendix C — `dnadiff`, Parsnp and wgMLST are different

> `dnadiff` compares two selected assemblies directly. Parsnp identifies a shared core across many genomes and supports a separate core-genome analysis. wgMLST counts allele differences using a named scheme. Results can be directionally compared, but the numbers cannot be converted into one another.

### Appendix D — Statistical methods and repeated observations

Many participants contribute repeated episodes. Formal episode-level comparisons should therefore account for within-participant dependence using an appropriate approach, such as a mixed model, clustered standard errors, generalised estimating equations or a participant-level analysis.

`14_genotype_phenotype_model.R` is the only modelling script named in the presentation because its exploratory association design is discussed. The sparse UTI outcome means that model convergence, singular fits, separation, multiple testing and overfitting require explicit reporting. Exploratory results must not be described as confirmatory.

### Appendix E — Claim-to-script traceability

Use this map to explain which code directly supports each method or result shown in the slides:

| Claim or displayed method | Script(s) that directly support it |
|---|---|
| Clinical cleaning and UTI/Not_UTI classification | `00a_load_clean_clinical.R`; `00b_classify_episodes.R` |
| Assembly reconciliation, implemented QC and canonical selection | `00_make_assembly_metadata.r`; `12a_wgs_qc.R` |
| VFDB-derived detected/not-detected matrix | `02_gene_presence_analysis.R` |
| Preferred provider/local-fallback MLST integration | `06_MLST.R` |
| Provenance-verified pairwise `dnadiff` and Mash comparison | `11_compare_strains.R` and its helper |
| Clinical, VF and MLST episode-level integration | `22_vf_build_analysis_dataset.R` |
| Mixed-canonical longitudinal reconstruction | `24_vf_longitudinal_dynamics.R` |
| Rebuilt Longcycler-only sensitivity timelines | `scripts/rebuild_longcycler_sensitivity.R` |
| Separate core-genome and pangenome context, where discussed | `12b_core_snp.R`; `12c_panaroo.R` |
| Exploratory association modelling and its cautions | `14_genotype_phenotype_model.R` |

Other scripts may be mentioned briefly when they clarify an upstream assumption or supporting diagnostic, but they are not presented as headline-result stages unless their outputs are actually used in the slide claim.

### Appendix F — Missing provenance

Do not invent missing details. Ask for or recover:

- sequencing platform, flow cell, library kit and run information;
- basecaller version, model and read-QC/coverage summaries;
- Longcycler and Flye versions, parameters and polishing steps;
- ABRicate and VFDB version/date;
- SeqSphere/MLST scheme version and provider manifest;
- MUMmer, Mash, Parsnp, `snp-dists`, Prokka and Panaroo versions;
- R package versions/session information; and
- citation or calibration for the operational SNP boundary.

## Rapid lecturer Q&A

**Why not simply delete Flye?**

> The executed pipeline selected Flye only as fallback. Excluding it is useful as a sensitivity check, but deleting it from the main record would conceal how the analysis was actually run and could introduce selection.

**Does one assembler make the analysis unbiased?**

> No. It standardises one technical choice. Sequencing quality, polishing, assembly quality, sampling and missing episodes can still affect results.

**Why use both `dnadiff` and Parsnp?**

> They answer related but distinct questions. Pairwise `dnadiff` supports direct within-participant assembly comparisons; Parsnp provides a shared-core population analysis across many genomes.

**Can the project prove persistence?**

> It can identify close genomic similarity between sampled timepoints. Continuous persistence between those timepoints remains an inference because intermediate bacterial populations were not continuously observed.

**Can the project prove what caused UTI?**

> Not from these genomic summaries alone. Causal inference would need a design and model addressing timing, measured confounders, repeated observations, treatment indication and reverse causality.

**How do I know the corrected SNP numbers are tied to the current assemblies?**

> The current pairwise output contains exactly 963 unordered within-participant comparisons from the 556 selected canonical genomes. Each comparison records both endpoint FASTA paths and SHA-256 hashes, and a cached report is accepted only when its provenance sidecar matches those current inputs. No report with a mismatched signature is reused.

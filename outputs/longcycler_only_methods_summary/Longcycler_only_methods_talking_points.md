# Longcycler-only methods slideshow: detailed talking points

Purpose: these notes are for explaining the Longcycler-only methods slideshow to a supervisor. The framing is deliberately conservative: it explains what the current analysis can support from available generated files, where the denominators come from, and why the primary defensible route uses one assembler rather than mixing Longcycler and Flye assemblies for strain-resolution claims.

## Opening framing

Suggested opening:

> This is a conservative methods framing. The existing mixed canonical analysis contains 556 genome-linked rows, but for the primary defensible route I restrict the genomic analysis to 532 QC-passing selected Longcycler rows from 161 participants. The reason is not that Flye is unusable, but that using one assembler avoids introducing assembler-driven technical heterogeneity when I am making strain, lineage, virulence-factor, and longitudinal comparisons.

Key points to establish early:

- The mixed canonical set is still useful context: it has 556 WGS/VF rows and reflects the broader generated pipeline output.
- The primary analysis denominator is narrower: 532 rows where the selected canonical assembly is Longcycler and QC-passing.
- The excluded 24 rows are Flye fallback rows. They are excluded for methodological consistency, not because they are necessarily wrong.
- This single-assembler decision is most important when discussing small genetic-distance thresholds, especially the strict same-strain rule of `≤25 SNPs`.
- Clinical labels are not being redefined because of assembler choice. The clinical `UTI` versus `Not_UTI` labels remain the clinical framework; the assembler restriction only changes which genome-linked rows are included in the primary genomic analysis.

Useful short version:

> I am separating clinical classification from genomic defensibility. Clinically, the labels stay the same. Genomically, I am using a single QC-passing assembler route so that any downstream differences are less likely to be artefacts of assembly choice.

## Important denominator language

Do not describe the counts as one simple shrinking dataset. They are different units:

- `rows` or `episodes` are sample-level observations.
- `participants` are residents/individuals.
- `genes` are VF feature columns, not samples.
- `STs` are lineage categories among typed rows.
- `transitions` are consecutive within-participant episode pairs.
- `strict same-strain transitions` are a subset of transitions supported by SNP evidence.

Safe wording:

> These counts are not all the same kind of denominator. The 532 rows are genome-linked episode rows; the 161 participants are people; the 227 VF genes are columns/features; and the 371 transitions are within-resident time-adjacent pairs rebuilt after applying the Longcycler-only restriction.

Why this matters:

- It prevents accidental comparison of incompatible units, such as saying 227 VF genes is a reduction from 532 rows.
- It makes clear that the longitudinal analysis has its own denominator because transitions are created from repeated observations within individuals.
- It explains why the same participant can contribute more than one row and more than one transition.

## Slide 1: Methods choice

Slide title: `Longcycler-only is the clean primary denominator`

Slide message:

This slide explains the methodological choice: the current generated analysis has a mixed canonical set, but the primary supervisor-facing route uses only selected, QC-passing Longcycler assemblies.

### What the title means

`Longcycler-only is the clean primary denominator` means:

- The primary genomic denominator is not every genome-linked row available.
- It is the subset where the selected canonical assembly was produced by Longcycler and passed QC.
- This makes the methods easier to defend because all included assemblies come from the same assembly route.

Suggested wording:

> I am using Longcycler-only as the clean primary denominator because it avoids making strain-resolution conclusions across a mixture of assembly methods. The 556-row mixed set remains useful background, but the primary claims should come from the 532 selected QC-passing Longcycler rows.

### Subtitle explanation

Slide subtitle:

`The analysis keeps one assembler for strain, lineage, and VF summaries; the 556-row mixed canonical set remains context.`

What this means:

- The mixed canonical set has 556 WGS/VF rows.
- That mixed set was Longcycler-preferred, but included 24 Flye fallback rows.
- The Longcycler-only analysis removes those Flye fallback rows.
- The mixed set can still be mentioned as context or sensitivity, but the cleaner primary denominator is 532 Longcycler rows.

Safe wording:

> I am not deleting the mixed analysis from the project. I am treating it as context, while using the Longcycler-only subset as the main defensible analysis route.

### Panel: “Why this helps”

Slide text:

`Single-assembler denominator reduces technical variation from assembly choice.`

Talking point:

Assemblers can differ in how they resolve repeats, structural complexity, polishing errors, and locus-level sequence calls. If a study is asking whether two isolates are closely related, or whether a small genetic threshold supports persistence, it is safer to avoid mixing assembler outputs unless assembler equivalence has been explicitly validated for that purpose.

Safe wording:

> A single assembler does not remove every technical source of variation, but it removes one avoidable source: differences introduced by using different assembly algorithms for different rows.

Slide text:

`It avoids defending equivalence between assemblers for small distance thresholds.`

Talking point:

This is especially important because the analysis uses a strict same-strain threshold of `≤25 SNPs`. When a threshold is small, even modest technical differences can matter. The point is not that Flye would necessarily change the result, but that defending mixed-assembler equivalence would require extra validation.

Safe wording:

> Because my strict same-strain threshold is small, I do not want the defence of the analysis to depend on assuming that Longcycler and Flye behave identically for all included isolates.

### Panel: “What changes”

Slide text:

`Primary genomic denominator becomes 532 Longcycler rows from 161 participants.`

Source and filter:

- Source file: `results/qc/canonical_assembly_selection.csv`
- Filter: selected canonical rows where `selected_canonical == TRUE`, `QC_PASS == TRUE`, and `assembler == "longcycler"` after case-normalisation.
- Interpretation: these are the genome-linked rows retained for the primary single-assembler analysis.

Talking point:

The primary denominator is therefore 532 episode-level rows, not 532 participants. These rows come from 161 participants, meaning some participants have repeated observations.

Slide text:

`The 24 Flye fallback rows are excluded rather than treated as interchangeable.`

Talking point:

The mixed canonical set has 556 rows. The Longcycler-only set has 532 rows. The difference is 24 rows. These 24 are Flye fallback rows that were included in the mixed canonical context but excluded from the primary single-assembler route.

Safe wording:

> Flye is not being labelled bad or invalid. I am simply choosing not to mix assembler routes in the primary analysis because that is easier to justify for strain-resolution work.

### Panel: “What does not change”

Slide text:

`Clinical UTI labels are unchanged and applied independently of assembler choice.`

Talking point:

The UTI label is a clinical classification. Restricting to Longcycler-only rows does not redefine UTI. It only determines which genome-linked episodes are available for genomic summaries.

Safe wording:

> The clinical definition is not being reverse-engineered from the genome data. The genome restriction is applied after clinical status has already been defined.

Slide text:

`ST remains lineage context; ≤25 SNPs remains the strict same-strain rule.`

Talking point:

MLST sequence type is useful for broad lineage context, for example ST131 or ST73. But ST is not high-resolution enough to prove that two isolates are the same strain. For strict persistence/same-strain claims, the analysis uses SNP evidence, with `≤25 SNPs` treated as strong same-strain support.

Safe wording:

> ST tells me whether isolates sit in the same broad lineage. SNP distance is what I use for the strict same-strain interpretation.

### Bottom framing note

Slide text:

`Discussion framing: conservative methods decision, not a claim that Flye assemblies are unusable.`

This sentence is important. It prevents the supervisor from hearing the slide as an attack on Flye. The methodological claim is narrower:

- Mixed assemblers are harder to defend for small-threshold genomic comparisons.
- Longcycler-only is cleaner and more conservative.
- Flye may still be useful in other contexts, but it is not part of the primary denominator here.

## Slide 2: Denominator ladder

Slide title: `The single-assembler analysis retains 532 genome-linked episodes`

Slide message:

This slide shows the main counts used in the Longcycler-only methods framing. The key thing to say is that each count has a specific unit, source file, filter, denominator, and interpretation.

Suggested opening for the slide:

> This slide is a denominator ladder, but not every number is the same kind of object. I am moving from the mixed canonical context into the Longcycler-only analysis, then into separate downstream layers: clinical status, VF genes, MLST, longitudinal transitions, and strict SNP-supported transitions.

### 01. Raw canonical set: 556

Visible text:

`Raw canonical set: 556`

`Longcycler-preferred + Flye fallback`

`current generated canonical analysis`

Source:

- `results/vf/vf_analysis_ready.csv`

Filter/definition:

- Current generated canonical WGS/VF-ready rows.
- This is the mixed canonical analysis set.
- It is Longcycler-preferred but includes Flye fallback rows where needed.

Denominator:

- 556 genome-linked episode rows.

Interpretation:

This is the broad current generated genomic context. It is not the primary single-assembler denominator, but it explains where the Longcycler-only set comes from.

Safe wording:

> The 556 rows are the current mixed canonical WGS/VF-ready context. This is the larger generated set, but it includes a small number of Flye fallback assemblies.

### 02. Longcycler-only set: 532

Visible text:

`Longcycler-only set: 532`

`QC-passing selected Longcycler rows`

`24 Flye fallback rows excluded`

`canonical_assembly_selection.csv`

Source:

- `results/qc/canonical_assembly_selection.csv`

Filter/definition:

- `selected_canonical == TRUE`
- `QC_PASS == TRUE`
- `assembler == "longcycler"` after normalising case.
- No failed, missing, non-selected, or rescued Longcycler assemblies are added back.

Denominator:

- 532 genome-linked episode rows.

Interpretation:

This is the primary defensible single-assembler denominator.

Safe wording:

> The 532 rows are the selected, QC-passing Longcycler rows. This is the main denominator I would defend in the methods, because every included genome follows the same assembly route.

### 03. Participants: 161

Visible text:

`Participants: 161`

`retained Longcycler-linked participants`

`vf_analysis_ready.csv`

Source:

- Participant IDs linked through the Longcycler-only rows.
- Cross-checked against `results/vf/vf_analysis_ready.csv` after applying the Longcycler-only row restriction.

Filter/definition:

- Unique participants represented among the 532 Longcycler-only rows.

Denominator:

- 161 participants.

Interpretation:

This tells us how many residents/individuals are represented, not how many samples exist. Because some residents have multiple timepoints, participant count is lower than row count.

Safe wording:

> These are 532 rows from 161 participants, so repeated sampling is present. That repeated structure is why transition analyses are built within participant over time.

### 04. Clinical status: 16 | 516

Visible text:

`Clinical status: 16 | 516`

`UTI | Not_UTI using primary status`

`status joined to Longcycler rows`

Source:

- Clinical status fields joined to the Longcycler-only genomic rows.

Filter/definition:

- Among the 532 Longcycler-only rows:
  - `UTI_Status == "UTI"` gives 16 rows.
  - `UTI_Status == "Not_UTI"` gives 516 rows.

Denominator:

- 532 Longcycler-only genome-linked episode rows.

Formula:

- `16 + 516 = 532`

Interpretation:

This is the clinical status distribution within the primary Longcycler-only genomic denominator.

Safe wording:

> Within the Longcycler-only set, there are 16 UTI rows and 516 Not_UTI rows. These are row-level clinical classifications, not participant-level counts.

### 05. VF feature space: 227

Visible text:

`VF feature space: 227`

`binary VF gene columns retained`

`vf_analysis_ready.csv`

Source:

- `results/vf/vf_analysis_ready.csv`

Filter/definition:

- VF gene columns retained in the analysis-ready matrix.
- The feature range is from `fyuA` through `sigA`.

Denominator:

- 227 VF gene columns/features.

Interpretation:

This is not a sample count. It is the number of virulence-factor features available for presence/absence analysis across the retained rows.

Safe wording:

> The 227 number is a feature count. It means there are 227 binary VF gene columns that can be summarised across the Longcycler-only rows.

### 06. MLST lineage: 514 | 80

Visible text:

`MLST lineage: 514 | 80`

`typed rows | distinct STs`

`lineage context only`

`mlst_provider_preferred.csv`

Source:

- `results/mlst/mlst_provider_preferred.csv`

Filter/definition:

- Longcycler-only selected rows joined to preferred MLST calls.
- Non-missing typed ST rows are counted.

Denominator:

- 514 typed rows.
- 80 distinct sequence types.

Interpretation:

MLST gives broad lineage context. It helps describe the population structure and lineages represented, but it does not by itself prove same-strain persistence.

Safe wording:

> MLST is used as lineage context. There are 514 typed Longcycler-only rows and 80 distinct STs. If two isolates share an ST, that is compatible with relatedness, but not enough to call them the same strain.

Why 514 is lower than 532:

- Some Longcycler-only genome-linked rows do not have a usable/non-missing ST call in the preferred MLST output.
- This is normal for a typed-subset denominator and should be reported separately.

### 07. Longitudinal pairs: 371

Visible text:

`Longitudinal pairs: 371`

`rebuilt in pipeline time order`

`from 139 participants`

`vf_longitudinal_transitions.csv`

Source:

- `results/vf/vf_longitudinal_transitions.csv`
- Longcycler-only reconstruction using the pipeline time-order fields.

Filter/definition:

- Restricted to `cohort == "all"` when using transition outputs, to avoid double-counting repeated cohort subsets.
- Longcycler-only rows are ordered within each participant.
- Consecutive pairs are rebuilt after excluding Flye fallback rows.

Denominator:

- 371 consecutive Longcycler-only transitions.
- These transitions come from 139 participants with at least two retained Longcycler-only timepoints.

Formula:

- For each participant, count retained Longcycler-only timepoints.
- If a participant has `n` retained timepoints, that participant contributes `n - 1` consecutive transitions.
- Sum across participants.

Interpretation:

This is the denominator for longitudinal transition analyses. It is not simply the number of participants or the number of rows.

Important explanation:

The transitions are rebuilt rather than just filtered from the mixed transition table because removing a Flye fallback row can change adjacency. For example, if a participant has:

`Longcycler timepoint 1 → Flye fallback timepoint 2 → Longcycler timepoint 3`

then filtering the mixed transition table would remove both original adjacent transitions and fail to create the new Longcycler-only adjacency:

`Longcycler timepoint 1 → Longcycler timepoint 3`

That is why the Longcycler-only transition table must be rebuilt in pipeline time order.

Safe wording:

> The 371 transitions are not obtained by naively filtering the mixed transition table. They are rebuilt after applying the Longcycler-only restriction, because once a Flye timepoint is removed the neighbouring Longcycler timepoints may become consecutive in the primary analysis.

### 08. Strict SNP evidence: 116

Visible text:

`Strict SNP evidence: 116`

`≤25 SNP among rebuilt pairs`

`7/9 Not_UTI -> UTI strict`

`pairwise_metrics.csv`

Source:

- `results/strain_compare/pairwise_metrics.csv`

Filter/definition:

- Longcycler-only rebuilt consecutive transitions.
- Pairwise SNP metrics joined to those transitions.
- Strict same-strain evidence defined as `TotalSNPs <= 25` or the equivalent current SNP context label `Strong same strain`.

Denominator:

- 116 strict same-strain transitions among the rebuilt Longcycler-only consecutive transitions.

Focused UTI-transition readout:

- There are 9 `Not_UTI -> UTI` transitions in the Longcycler-only rebuilt transition set.
- Of these 9, 7 have strict same-strain SNP evidence.

Interpretation:

The 116 count is the strict SNP-supported same-strain subset. The `7/9` readout is a focused clinical transition result, showing that most Longcycler-only Not_UTI-to-UTI transitions with this definition are consistent with same-strain continuity under the strict SNP rule.

Safe wording:

> I am using SNP distance for the strict same-strain call. In the Longcycler-only longitudinal rebuild, 116 transitions meet the ≤25 SNP criterion. For the focused Not_UTI-to-UTI transition, 7 of 9 are strict same-strain by that SNP threshold.

Caution:

Do not overstate this as proof that colonisation caused UTI. It supports strain continuity across some clinical transitions, not causality.

## Slide 3: Methods boundary

Slide title: `Longcycler-only supports reproducible methods`

Slide message:

This slide defines what the Longcycler-only analysis supports and what it does not attempt to support. It separates factual reproducibility from biological interpretation and from causal claims.

Suggested opening:

> This slide is the boundary slide. It tells the supervisor what I am comfortable claiming from the Longcycler-only generated outputs, and what I am deliberately not recalculating from proxies.

### Subtitle explanation

Slide subtitle:

`Boundary: reproducibility choices are separate from biological interpretation.`

Talking point:

The single-assembler choice is a reproducibility and methods decision. It does not automatically create biological conclusions. The biological interpretation still needs careful wording and should not imply causality unless the design supports it.

Safe wording:

> This slide separates the technical denominator from the biological interpretation. The Longcycler-only restriction makes the genomic workflow cleaner, but it does not by itself prove mechanisms or causality.

### Panel: “Supported”

Slide text:

`VF presence/absence, MLST lineage context, core SNP strain context, and longitudinal transitions.`

Talking point:

The Longcycler-only set supports:

- VF presence/absence summaries across retained genome-linked rows.
- MLST lineage summaries among typed rows.
- SNP-based strain context for linked pairs where pairwise SNP metrics are available.
- Longitudinal transition summaries rebuilt within participant in pipeline time order.

Safe wording:

> The Longcycler-only set lets me summarise virulence-factor carriage, describe MLST lineages, evaluate SNP-supported strain continuity, and examine within-resident transitions over time.

Slide text:

`Strict persistence is reported with SNP evidence, not ST agreement alone.`

Talking point:

This prevents overcalling persistence. Two isolates with the same ST might still be different strains within a common lineage. Therefore, strict persistence/same-strain classification is based on SNP evidence.

Safe wording:

> I use ST to orient the lineage, but I do not use ST alone as proof of persistence. Strict persistence requires the SNP threshold.

### Panel: “Not recalculated”

Slide text:

`No antibiotic, demographic, host-factor, or wgMLST allele claims are inferred from proxies.`

Talking point:

This is a strong reproducibility boundary. If the current analysis does not contain the necessary antibiotic exposure data, host-factor data, demographic extract, or wgMLST allele-distance table, those claims should not be recomputed using substitute variables.

Safe wording:

> I am not trying to recreate antibiotic, demographic, host-factor, or wgMLST claims unless the required source data are available. I would rather mark those as not checkable than calculate them from proxies.

Slide text:

`The summary does not treat SNP and wgMLST distances as interchangeable.`

Talking point:

SNP distances and wgMLST allele distances are both genomic relatedness measures, but they are not the same measurement. They use different units, different calling processes, and potentially different error structures. They can be directionally compared, but one should not be substituted for the other.

Safe wording:

> SNP and wgMLST can both speak to relatedness, but a SNP threshold is not automatically equivalent to a wgMLST allele threshold.

### Panel: “Supervisor ask”

Slide text:

`Use Longcycler-only as the primary defensible methods denominator.`

Talking point:

This is the practical request. The proposed primary methods route is:

- Use 532 selected QC-passing Longcycler rows for primary genomic summaries.
- Use the 556-row mixed set as context or sensitivity only.
- Keep labels and downstream definitions explicit.

Safe wording:

> My suggested primary denominator is the Longcycler-only set. It is a small reduction in rows but a clearer methodological position.

Slide text:

`Keep the 556-row mixed canonical analysis only as context or sensitivity, if needed.`

Talking point:

This preserves the broader analysis without making it the main defensible result. If the supervisor wants, the mixed set can be shown as a sensitivity analysis to demonstrate whether headline patterns are stable when Flye fallback rows are included.

Safe wording:

> The mixed set does not need to disappear. It can remain as a sensitivity or context analysis, but I would not lead with it for strain-resolution claims.

### Bottom suggested wording

Slide text:

`Suggested wording: primary genomic analyses were restricted to QC-passing Longcycler assemblies to minimise assembler-induced technical heterogeneity.`

This is the most manuscript-like sentence on the slide. It is safe because it:

- says exactly what was done;
- gives a methodological reason;
- does not claim that Longcycler is universally best;
- does not claim that Flye is invalid;
- avoids causal or biological overreach.

## Slide 4: Analysis flow

Slide title: `How the results flow through the Longcycler-only analysis`

Slide message:

This slide shows the pathway from clinical episodes and assembly/QC records into the mixed canonical set, then into the Longcycler-only primary denominator, and finally into downstream analyses.

Suggested opening:

> This flowchart shows where each denominator enters the analysis. The clinical and assembly/QC inputs feed the mixed canonical context, then I apply the Longcycler-only restriction and rebuild the downstream summaries from that primary denominator.

### Box: Clinical status

Visible text:

`Clinical status`

`583 episodes`

`166 participants`

`18 UTI | 565 Not_UTI`

Source/definition:

- Clinical episode-level dataset used to define UTI versus Not_UTI status.

Denominator:

- 583 clinical episodes.
- 166 participants.
- 18 UTI episodes.
- 565 Not_UTI episodes.

Interpretation:

This is the broader clinical starting point before restricting to WGS/VF-ready Longcycler-only genomic rows.

Safe wording:

> The clinical denominator starts larger than the genomic denominator. There are 583 clinical episodes from 166 participants, with 18 UTI and 565 Not_UTI episodes before genome-linked assembly restrictions are applied.

### Box: Assembly and QC

Visible text:

`Assembly and QC`

`1,291 records`

`Longcycler/Flye alternatives`

`QC and canonical selection`

Source/definition:

- `results/qc/canonical_assembly_selection.csv`

Denominator:

- 1,291 assembly records.

Interpretation:

This is not the number of final samples. It is the assembly/QC universe, where there can be alternative assemblies for a sample, including Longcycler and Flye candidates. QC and canonical-selection logic decide which assembly is selected for downstream analysis.

Safe wording:

> The 1,291 records are assembly records, not final analytic rows. This reflects that samples can have alternative assembly records, and the pipeline then selects a canonical assembly after QC.

### Arrow into mixed canonical context

Meaning:

Clinical status and assembly/QC information are brought together to create the current genome-linked canonical analysis context.

Safe wording:

> The clinical labels and the selected assembly records come together in the WGS/VF-ready canonical analysis table.

### Box: Mixed canonical context

Visible text:

`Mixed canonical context`

`556 WGS/VF rows`

`162 participants`

`Longcycler-preferred + 24 Flye fallback`

Source/definition:

- `results/vf/vf_analysis_ready.csv`

Denominator:

- 556 WGS/VF-ready rows.
- 162 participants.
- Includes 24 Flye fallback rows.

Interpretation:

This is the current generated canonical set. It is broader than the primary Longcycler-only set and remains useful context.

Safe wording:

> The mixed canonical context contains 556 WGS/VF-ready rows from 162 participants. It is Longcycler-preferred, but it includes 24 Flye fallback rows.

### Arrow into primary Longcycler denominator

Meaning:

This arrow represents the key methodological restriction: from the mixed canonical context, retain only selected QC-passing Longcycler assemblies.

Safe wording:

> This is the main methods decision: I move from the mixed canonical context into the single-assembler Longcycler-only primary denominator.

### Box: Primary Longcycler denominator

Visible text:

`Primary Longcycler denominator`

`532 rows | 161 participants`

`16 UTI | 516 Not_UTI`

`selected canonical QC-passing Longcycler`

Source/definition:

- `results/qc/canonical_assembly_selection.csv`
- joined/checked against the analysis-ready VF and clinical status outputs.

Filter:

- selected canonical row;
- QC-passing;
- assembler is Longcycler.

Denominator:

- 532 genome-linked rows.
- 161 participants.
- 16 UTI rows.
- 516 Not_UTI rows.

Interpretation:

This is the primary denominator for supervisor-facing methods and downstream genomic summaries.

Safe wording:

> This is the denominator I would lead with: 532 selected QC-passing Longcycler rows from 161 participants, containing 16 UTI and 516 Not_UTI rows.

### Box: Analysis layers

Visible text:

`Analysis layers`

`MLST: 514 | 80 STs`

`VF matrix: 227 genes`

`SNP context: ≤25 = strict same strain`

Sources:

- MLST: `results/mlst/mlst_provider_preferred.csv`
- VF matrix: `results/vf/vf_analysis_ready.csv`
- SNP context: `results/strain_compare/pairwise_metrics.csv`

Interpretation:

This box separates the downstream genomic layers:

- MLST tells us lineage context.
- VF matrix tells us virulence-factor presence/absence.
- SNP context tells us high-resolution strain relatedness for pairwise comparisons.

Safe wording:

> Once I have the Longcycler-only denominator, I analyse it through separate layers: MLST for lineage, VF genes for gene carriage, and SNPs for strict strain-continuity evidence.

Important caution:

The three layers do not have identical denominators:

- 532 rows in the Longcycler-only denominator.
- 514 typed MLST rows.
- 227 VF gene columns.
- SNP evidence applies to eligible pairwise transitions, not every row.

### Box: Longitudinal rebuild

Visible text:

`Longitudinal rebuild`

`371 transitions`

`139 participants`

`rebuilt in pipeline time order`

Source:

- `results/vf/vf_longitudinal_transitions.csv`
- rebuilt after applying Longcycler-only eligibility.

Definition:

- Within each participant, retained Longcycler-only timepoints are ordered by pipeline time order.
- Consecutive transitions are rebuilt.
- Participants with only one retained timepoint do not contribute transitions.

Denominator:

- 371 transitions.
- 139 participants.

Interpretation:

This is the transition denominator for longitudinal analyses after applying the single-assembler restriction.

Safe wording:

> The longitudinal transition denominator is rebuilt after excluding Flye fallback rows. That gives 371 consecutive Longcycler-only transitions from 139 participants.

### Box: Focused readout

Visible text:

`Focused readout`

`9 Not_UTI→UTI`

`7 strict same-strain`

`among 116 strict pairs`

Sources:

- Clinical transition labels from the rebuilt Longcycler-only transition table.
- SNP pairwise metrics from `results/strain_compare/pairwise_metrics.csv`.

Definition:

- Count transitions where the previous row is `Not_UTI` and the next row is `UTI`.
- Among those, count transitions meeting strict same-strain SNP evidence.

Denominator:

- 9 Not_UTI-to-UTI transitions.
- 7 of those 9 are strict same-strain.
- 116 total strict same-strain transitions across the rebuilt Longcycler-only transition set.

Interpretation:

This is a focused clinical-transition result. It supports the idea that many observed Not_UTI-to-UTI transitions in this subset involve close SNP continuity, but it does not prove that carriage caused UTI.

Safe wording:

> In the Longcycler-only rebuild, there are 9 Not_UTI-to-UTI transitions, and 7 meet the strict ≤25 SNP same-strain criterion. I would phrase this as strain continuity across transition, not as causation.

### Boundary note

Visible text:

`Boundary: 24 Flye fallback rows excluded; antibiotic, wgMLST, demographic and causal claims are not recalculated here.`

Why it matters:

This prevents the flowchart from being interpreted as a full reproduction of every possible thesis claim. It is a methods and genomic-output flowchart, not an antibiotic, host-factor, wgMLST, or causal model reproduction.

Safe wording:

> This workflow is intentionally bounded. It does not recalculate antibiotic, demographic, host-factor, wgMLST, or causal claims from proxy variables.

## Likely supervisor questions and strong answers

### 1. “Why not use Flye?”

Answer:

> I am not saying Flye is unusable. The issue is defensibility and consistency. The mixed canonical set includes 24 Flye fallback rows, but for the primary strain-resolution analysis I would rather use one assembler throughout. That avoids having to defend assembler equivalence when using small genetic thresholds such as ≤25 SNPs. Flye can remain in the mixed set as context or sensitivity, but I would not lead with a mixed-assembler denominator.

Short version:

> It is a conservative consistency decision, not a rejection of Flye.

### 2. “Are you losing too much data by excluding Flye?”

Answer:

> The reduction is modest: from 556 mixed canonical rows to 532 Longcycler-only rows, so 24 rows are excluded. The participant count changes from 162 to 161. In exchange, the primary analysis becomes much cleaner methodologically because it avoids mixing assembler routes. For supervisor discussion, I think that trade-off is reasonable.

Useful numbers:

- Rows lost: `556 - 532 = 24`
- Participant count change: `162` to `161`
- Longcycler-only rows retained: `532 / 556`, approximately 95.7% of the mixed canonical rows.

Safe wording:

> I retain the large majority of the generated canonical rows while removing a technical complication.

### 3. “Why does the UTI count change from 18 to 17 to 16?”

Answer:

> These are different denominators. The 18 UTI episodes are in the broader clinical-status denominator of 583 episodes. When the analysis is restricted to genome-linked canonical outputs, the UTI count can change because not every clinical episode necessarily has an eligible genome-linked row. Under the Longcycler-only primary denominator, the retained rows include 16 UTI and 516 Not_UTI rows.

Important nuance:

- The slide deck explicitly shows `18 UTI | 565 Not_UTI` at the clinical-starting level.
- The Longcycler-only denominator shows `16 UTI | 516 Not_UTI`.
- If a 17-row UTI count appears in another current output, describe it as an intermediate or alternate genomic denominator and verify its exact source before presenting it as final.

Safe wording:

> The UTI count changes because the denominator changes, not because the clinical definition is being rewritten.

### 4. “Why is ST not enough to call same strain?”

Answer:

> ST is useful but low-resolution. Two isolates can share the same ST and still be different strains within that lineage. That is especially true for common lineages. Therefore, I use ST for lineage context and SNP distance for strict same-strain evidence.

Safe wording:

> Same ST is compatible with relatedness; it is not proof of same-strain persistence.

### 5. “What does ≤25 SNPs mean?”

Answer:

> In this analysis, ≤25 SNPs is used as the strict same-strain threshold for consecutive linked pairs. If a rebuilt Longcycler-only transition has a pairwise SNP distance of 25 or fewer SNPs, I classify it as strong same-strain evidence. This is a conservative high-resolution rule compared with relying on ST alone.

Caution:

> The threshold supports strain continuity; it does not by itself prove direction, transmission, or causation.

### 6. “Why are wgMLST and SNP distances not interchangeable?”

Answer:

> They measure relatedness differently. SNP distances count nucleotide differences in a defined comparison framework, while wgMLST distances count allele differences across a whole-genome MLST scheme. The units, calling process, missingness, and error structures differ. They can be directionally corroborative, but a SNP threshold cannot be assumed to equal a wgMLST allele threshold.

Safe wording:

> I would use wgMLST only if the actual allele-distance table and scheme details are available. I would not convert SNP results into wgMLST claims.

### 7. “Are you making causal claims?”

Answer:

> No. The Longcycler-only analysis supports descriptive and reproducibility claims: which isolates are present, which lineages and VF genes are observed, and whether consecutive pairs show strict SNP-supported strain continuity. It does not establish that antibiotics, colonisation, or any host factor caused UTI or replacement unless the design and model specifically support that.

Safe wording:

> I am careful to say “consistent with” or “supports strain continuity,” not “causes.”

### 8. “How would this be described in the methods section?”

Possible methods wording:

> Primary genomic analyses were restricted to selected, QC-passing Longcycler assemblies to minimise assembler-induced technical heterogeneity. The broader mixed canonical dataset contained 556 WGS/VF-ready rows and was Longcycler-preferred with Flye fallback assemblies. For the primary single-assembler analysis, Flye fallback rows were excluded, retaining 532 Longcycler-derived genome-linked rows from 161 participants. Clinical UTI status was applied independently of assembler choice. MLST was used for lineage context, VF genes were analysed as binary presence/absence features, and strict same-strain continuity across consecutive within-participant transitions was defined using pairwise core-genome SNP evidence with a threshold of ≤25 SNPs.

If the supervisor wants a shorter version:

> To avoid mixed-assembler effects, the primary genomic analysis used only selected QC-passing Longcycler assemblies. This retained 532 rows from 161 participants and was used for VF, MLST, SNP-context, and longitudinal transition summaries.

## Safe wording table

| Topic | Safer wording | Avoid saying |
|---|---|---|
| Flye | “Flye fallback rows were excluded from the primary analysis to keep a single assembler.” | “Flye is wrong” or “Flye cannot be used.” |
| Longcycler-only | “Primary defensible methods denominator.” | “The only true dataset.” |
| Mixed canonical set | “Context or sensitivity analysis.” | “Invalid analysis.” |
| ST | “Lineage context.” | “Same ST proves same strain.” |
| SNP threshold | “≤25 SNPs gives strict same-strain evidence.” | “≤25 SNPs proves causation.” |
| Not_UTI→UTI | “7 of 9 transitions show strict same-strain continuity.” | “Colonisation caused UTI in 7 cases.” |
| wgMLST | “Not recalculated without allele-distance table and scheme details.” | “SNPs and wgMLST are equivalent.” |
| Antibiotics/host factors | “Not recalculated from proxies.” | “Antibiotics caused replacement/UTI” without proper model/data. |

## Factual reproducibility versus interpretation versus limitation

### Factual reproducibility claims supported here

- The mixed canonical WGS/VF-ready context contains 556 rows.
- The Longcycler-only primary denominator contains 532 selected QC-passing Longcycler rows.
- These 532 rows come from 161 participants.
- Within those rows, there are 16 UTI and 516 Not_UTI rows.
- The VF matrix contains 227 binary VF gene columns.
- The MLST typed subset contains 514 typed rows and 80 distinct STs.
- The Longcycler-only longitudinal rebuild contains 371 transitions from 139 participants.
- The strict SNP-supported transition subset contains 116 transitions.
- The focused Not_UTI-to-UTI readout contains 9 transitions, 7 of which are strict same-strain by the SNP threshold.

### Scientific interpretation that is reasonable but should be cautious

- The Longcycler-only restriction improves methodological consistency for genomic comparisons.
- ST distributions can describe the lineage structure of the retained isolates.
- VF presence/absence can be summarised across retained rows.
- SNP-supported continuity can be used to describe persistence-like patterns across consecutive timepoints.
- The `7/9` Not_UTI-to-UTI result is consistent with strain continuity across many of those transitions.

### Limitations to state clearly

- Excluding Flye fallback rows slightly reduces the dataset.
- Longcycler-only results may differ from the mixed canonical analysis if excluded rows are not random.
- MLST is not sufficient for strict same-strain inference.
- SNP and wgMLST distances should not be treated as interchangeable.
- Antibiotic, demographic, host-factor, wgMLST allele-distance, and causal claims are not recalculated here because they require specific source data and model details.
- Repeated observations within participants mean analyses should respect within-person dependence where formal modelling is performed.

## Final 60-second summary

If you need to explain the whole slideshow quickly:

> The main methodological decision is to use Longcycler-only as the primary genomic denominator. The current generated mixed canonical set has 556 WGS/VF rows, but it includes 24 Flye fallback rows. To avoid mixing assemblers when making strain-resolution claims, I restrict the primary analysis to 532 selected QC-passing Longcycler rows from 161 participants. Within that set there are 16 UTI and 516 Not_UTI rows. MLST is used for lineage context, with 514 typed rows and 80 STs, and the VF matrix contains 227 binary gene features. For longitudinal analysis, I rebuild transitions after applying the Longcycler-only restriction, giving 371 consecutive transitions from 139 participants. Strict same-strain evidence is based on ≤25 SNPs, giving 116 strict transitions overall; among 9 Not_UTI-to-UTI transitions, 7 meet that strict SNP criterion. I would present this as a conservative, reproducible methods framework, not as a causal claim and not as a statement that Flye is unusable.


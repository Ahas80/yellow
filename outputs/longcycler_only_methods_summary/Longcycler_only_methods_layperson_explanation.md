# Longcycler-only slideshow: simple explanation ladder

Purpose: this is the plain-English version of the Longcycler-only slideshow explanation. Use it when explaining the work to someone who is not deep in genomics, bioinformatics, or statistics. It starts very simple, then gradually adds detail so you can decide how far to go depending on the audience.

## The one-sentence explanation

> I am using one consistent genome-assembly method, Longcycler, for my main analysis so that any differences I see between samples are more likely to be biological rather than caused by mixing different computer methods.

## The 30-second explanation

> My full generated dataset has 556 genome-linked rows. Most were assembled with Longcycler, but 24 used Flye as a fallback. For my main analysis, I am being conservative and only using the 532 rows that passed QC and were assembled with Longcycler. This keeps the method consistent. I am not saying Flye is bad; I am just avoiding a mixed-method dataset when making strain-level comparisons. This gives me 532 rows from 161 participants, including 16 UTI rows and 516 Not_UTI rows.

## The simple analogy

Imagine you are comparing photographs to decide whether two pictures show the same person over time.

- If all photographs were taken using the same camera and same settings, the comparison is easier to defend.
- If most were taken with one camera but a few were taken with another camera, the comparison might still be okay, but someone could ask whether differences are due to the person or the camera.
- My decision is like saying: for the main analysis, I will use only the photos taken with the same camera.

In this analogy:

- The “camera” is the genome assembler.
- Longcycler and Flye are two different “cameras”.
- The “photos” are bacterial genome assemblies.
- The question is whether bacteria at different timepoints are closely related or the same strain.

Safe phrase:

> I am not saying the second camera is broken. I am saying that for the cleanest comparison, I want one camera in the main analysis.

## The explanation ladder

### Level 1: What is the big decision?

The big decision is:

> Use only Longcycler assemblies for the main analysis.

Why?

> Because it keeps the analysis consistent.

What changes?

> The dataset goes from 556 mixed genome-linked rows to 532 Longcycler-only rows.

What does not change?

> The clinical labels do not change. A UTI row is still a UTI row. A Not_UTI row is still a Not_UTI row.

### Level 2: Why does using one assembler matter?

A genome assembler is a computer method that takes DNA sequencing reads and builds them into a more complete genome sequence.

Different assemblers can make slightly different decisions. That matters when you are asking a sensitive question like:

> Are these two bacterial samples basically the same strain?

If you mix assemblers, a supervisor could reasonably ask:

> Are the differences real bacterial differences, or are they partly caused by using different assembly methods?

Your answer is:

> To avoid that concern, my primary analysis uses one assembler only: Longcycler.

### Level 3: What are the main numbers?

Use these numbers in this order:

| Number | Plain-English meaning |
|---:|---|
| 583 | Clinical episodes before genomic filtering. Think “clinical rows/visits”. |
| 166 | Participants in the broader clinical dataset. Think “people”. |
| 18 / 565 | Clinical UTI / Not_UTI episodes before restricting to genome-linked Longcycler rows. |
| 556 | Genome-linked rows in the mixed canonical analysis. Mostly Longcycler, with 24 Flye fallback rows. |
| 532 | Main Longcycler-only genome-linked rows. This is the clean primary denominator. |
| 161 | Participants represented in the 532 Longcycler-only rows. |
| 16 / 516 | UTI / Not_UTI rows within the Longcycler-only dataset. |
| 227 | Virulence-factor gene features. These are columns, not samples. |
| 514 / 80 | 514 rows have MLST typing; these contain 80 different STs. |
| 371 | Consecutive within-person transitions after rebuilding the timeline using Longcycler-only rows. |
| 116 | Transitions with strict SNP evidence for same strain. |
| 9 | Not_UTI-to-UTI transitions in the Longcycler-only rebuild. |
| 7 | Of those 9 Not_UTI-to-UTI transitions, 7 are strict same-strain by SNP evidence. |

The most important sentence:

> These numbers are different kinds of things: rows, people, genes, lineages, and transitions. They should not all be treated as one shrinking denominator.

## How to explain the dataset without confusing people

Use this:

> I start with clinical episodes, then I look at which of those have usable genome data, then I restrict the main genome analysis to Longcycler-only rows.

Then explain:

- Clinical episodes are the medical/clinical observations.
- Genome-linked rows are the clinical observations that also have usable bacterial genome data.
- Longcycler-only rows are the genome-linked rows assembled with the same method.
- Transitions are before-and-after comparisons within the same participant.

Very simple version:

> First: who had what clinical status?  
> Second: which episodes had usable genome data?  
> Third: which genome rows were assembled using Longcycler?  
> Fourth: among repeated samples from the same person, do the bacteria look like the same strain over time?

## Slide-by-slide layperson explanation

## Slide 1: Why Longcycler-only?

Main message:

> This slide explains why I chose one consistent genome-assembly method for the main analysis.

Say:

> The larger analysis has 556 genome-linked rows, but 24 of these used Flye as a fallback. To keep the main analysis cleaner, I use the 532 rows assembled with Longcycler. This helps avoid the criticism that small genetic differences might be due to mixing assembly methods.

If someone asks, “Is Flye bad?”

Say:

> No. I am not saying Flye is bad. I am saying that for my primary analysis, it is cleaner to avoid mixing methods.

If someone asks, “Why does this matter?”

Say:

> Because I am making strain-level comparisons. When you are asking whether bacteria are very closely related, even small technical differences can matter.

## Slide 2: What do all the numbers mean?

Main message:

> This slide shows the main denominators, but the numbers are not all measuring the same thing.

Simple explanation:

> Some numbers are sample rows, some are people, some are gene features, and some are before-and-after comparisons. So I should explain each number with its unit.

Go through the numbers like this:

### 556 mixed rows

> This is the larger genome-linked dataset. It includes mostly Longcycler rows plus 24 Flye fallback rows.

### 532 Longcycler-only rows

> This is my main clean dataset. These are selected, QC-passing rows assembled with Longcycler.

### 161 participants

> These 532 rows come from 161 people. Some people have more than one sample.

### 16 UTI and 516 Not_UTI rows

> Within the Longcycler-only dataset, 16 rows are UTI and 516 are Not_UTI.

Important clarification:

> These are row-level counts, not person-level counts.

### 227 VF genes

> These are virulence-factor gene features. They are like a checklist of 227 possible genes.

Analogy:

> If each bacterial sample is a row in a spreadsheet, the VF genes are columns saying present or absent.

### 514 typed rows and 80 STs

> MLST gives a broad bacterial family or lineage label. In the Longcycler-only set, 514 rows have usable MLST typing, and those rows contain 80 different sequence types.

Simple caution:

> ST is useful for broad grouping, but it is not precise enough by itself to prove two samples are the same strain.

### 371 transitions

> These are before-and-after comparisons within the same person.

Example:

> If one participant has three retained timepoints, they contribute two transitions: timepoint 1 to 2, and timepoint 2 to 3.

Why they are rebuilt:

> When I remove Flye fallback rows, I need to rebuild the timeline using only Longcycler rows. I cannot just filter the old transition table, because removing a middle timepoint can create a new Longcycler-to-Longcycler comparison.

### 116 strict same-strain transitions

> These are transitions where the bacteria are very close genetically by the SNP rule.

Simple version:

> Out of the Longcycler-only transitions, 116 have strong evidence that the strain stayed the same.

### 9 Not_UTI-to-UTI transitions, 7 strict same-strain

> There are 9 cases where a participant moves from Not_UTI to UTI in consecutive Longcycler-only timepoints. In 7 of those 9, the bacteria are close enough by SNPs to count as strict same-strain.

Safe wording:

> This suggests strain continuity across many of those transitions, but it does not prove that colonisation caused the UTI.

## Slide 3: What can and cannot be claimed?

Main message:

> This slide protects the analysis from overclaiming.

Say:

> The Longcycler-only dataset can support VF summaries, MLST lineage summaries, SNP strain-context summaries, and longitudinal transition summaries. But I am not using it to recalculate antibiotic effects, demographics, host factors, wgMLST allele distances, or causal claims.

Simpler version:

> This analysis is good for describing the bacteria and their relatedness over time. It is not enough by itself to prove what caused UTI or replacement.

Use this if asked about causality:

> I can say “consistent with strain continuity” or “supports same-strain persistence.” I should not say “this caused UTI” unless I have the right causal model and data.

Use this if asked about wgMLST:

> SNP distances and wgMLST allele distances are related ideas, but they are not the same measurement. I should not treat them as interchangeable.

## Slide 4: How does everything flow?

Main message:

> This slide shows the analysis pipeline from clinical data to genome data to Longcycler-only results.

Explain it slowly:

### Step 1: Clinical status

> I begin with 583 clinical episodes from 166 participants. These include 18 UTI and 565 Not_UTI episodes.

### Step 2: Assembly and QC

> Separately, there are 1,291 assembly/QC records. These are not final samples; they are assembly records, because some samples can have more than one assembly attempt or assembler option.

### Step 3: Mixed canonical context

> After choosing canonical assemblies and linking them to the analysis-ready genome data, the mixed context has 556 WGS/VF rows from 162 participants. This includes 24 Flye fallback rows.

### Step 4: Primary Longcycler denominator

> For the main analysis, I keep only the selected QC-passing Longcycler rows. That gives 532 rows from 161 participants, with 16 UTI and 516 Not_UTI rows.

### Step 5: Analysis layers

> From those 532 rows, I look at different analysis layers: MLST for lineage, VF genes for virulence-factor carriage, and SNP distances for strict same-strain evidence.

### Step 6: Longitudinal rebuild

> Then I rebuild the timeline within each participant using only Longcycler rows. This gives 371 before-and-after transitions from 139 participants.

### Step 7: Focused Not_UTI-to-UTI readout

> Finally, I look specifically at Not_UTI-to-UTI transitions. There are 9 of these, and 7 have strict same-strain evidence.

Final boundary sentence:

> This flowchart is about the genomic analysis route. It does not recalculate antibiotics, demographics, host factors, wgMLST allele distances, or causal effects.

## The “perfect explanation” script

You can say this almost word-for-word:

> The main thing I’m trying to do here is make the analysis easier to defend. The full generated genome-linked dataset has 556 rows, but it is a mixed canonical set because 24 rows use Flye fallback assemblies. For the primary analysis, I restrict to 532 selected, QC-passing Longcycler rows from 161 participants. That way, when I compare bacterial strains over time, I am not also mixing assembly methods.
>
> The clinical labels stay the same. In the Longcycler-only set, there are 16 UTI rows and 516 Not_UTI rows. I then use this same clean denominator for the main genomic summaries: 227 virulence-factor genes, 514 MLST-typed rows across 80 STs, and longitudinal strain comparisons.
>
> For the longitudinal part, I rebuild the timeline after excluding Flye fallback rows. That gives 371 consecutive transitions from 139 participants. Of those, 116 meet the strict same-strain SNP threshold of 25 SNPs or fewer. When I focus specifically on Not_UTI-to-UTI transitions, there are 9, and 7 of those show strict same-strain evidence.
>
> The cautious interpretation is that this supports strain continuity in many of these transitions. I am not claiming causation, and I am not recalculating antibiotic, demographic, host-factor, or wgMLST claims without the proper source data.

## Phrases to use

Use these:

- “Conservative methods choice.”
- “Primary defensible denominator.”
- “Single-assembler analysis.”
- “Avoids assembler-driven technical variation.”
- “Flye is excluded from the primary route, not labelled unusable.”
- “ST is lineage context.”
- “SNP distance gives stricter strain-level evidence.”
- “Consistent with strain continuity.”
- “Not a causal claim.”
- “Not recalculated from proxy data.”

## Phrases to avoid

Avoid these:

- “Flye is wrong.”
- “Longcycler is always the best.”
- “ST proves it is the same strain.”
- “SNPs and wgMLST are the same.”
- “Colonisation caused the UTI.”
- “Antibiotics caused replacement.”
- “These results prove causality.”
- “All the numbers are the same denominator.”

## If you get flustered, return to this

> I am using Longcycler-only as my main route because it is cleaner and easier to defend. It keeps one assembly method for the main genomic comparisons. The trade-off is losing 24 Flye fallback rows, but I still retain 532 of 556 mixed canonical rows, which is most of the data. The analysis can support descriptive genomic and strain-continuity findings, but I am careful not to overclaim causality or recalculate results that need data I do not have.


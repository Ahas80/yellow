# Detailed Presenter Script

**Deck:** Virulence-factor profiling of longitudinal urinary *E. coli* isolates

**Framing:** VF analysis of longitudinal urinary *E. coli* isolates first; UTI versus `Not_UTI` only as an exploratory clinical annotation.

## Avoid Saying
- This is a recurrent UTI presentation.
- ASB versus UTI is the current primary comparison.
- lpf is significantly associated with UTI in the current final analysis.
- Stable VF profile proves host-state mechanism.
- No VF gene change means nothing bacterial changed.
- AMR is a main global association result in this deck.
- Same sequence type proves same strain.

## Spoken Slides

### Slide 1. Virulence-factor profiling of longitudinal urinary E. coli isolates (1:00)
**Primary message:** Reset the frame: this is a VF-pipeline and longitudinal-isolate review, not a recurrent-UTI-only talk.
**Visual walkthrough:** Point to the three question cards, then the four denominator tiles at the bottom.

**Speaker script:**

I want to frame this review around the current virulence-factor pipeline for longitudinal urinary E. coli isolates. The clinical labels matter, but the main object of study is the VF repertoire measured from repeated WGS-linked isolates.

The three questions structure the talk. First: what VF repertoire is present in the current VF-ready dataset? Second: how stable are those profiles over time within the same resident? Third: once the pipeline is established, how can we overlay UTI versus Not_UTI status without overstating sparse clinical associations?

The bottom numbers are the key anchors for the rest of the review: 556 VF-ready episodes, 227 VF gene columns, 394 longitudinal comparisons, and 17 UTI annotations in the VF/model-ready clinical overlay.

**Emphasise:**
- Say VF analysis first; clinical status overlay second.
- Make clear that the deck is about the current pipeline and what it can safely support.

**Avoid:** Do not introduce this as a recurrent UTI findings deck or as a disease-association result.

**Transition:** I’ll start with why repeated-isolate VF profiling gives us a different question from a single clinical-status comparison.

### Slide 2. The YELLOW routine in one plain-English loop (1:00)
**Primary message:** Give a newcomer the whole behind-the-scenes loop before technical results appear.
**Visual walkthrough:** Move left to right across the six boxes, then use the three lower boxes to separate measurement, comparison, and clinical overlay.

**Speaker script:**

This slide is the simple mental model for the whole project. We start with a clinical episode: a urine-sampling moment with symptoms, culture context, and later clinical status.

From that episode, a urinary E. coli isolate is selected and sequenced. The sequencing data are assembled into a genome representation, and that genome is screened for virulence-factor genes.

The key data object is then a VF gene row: for each isolate, each virulence-factor gene is marked detected or not detected. Those rows are summarised into profile outputs such as modules, burden, Jaccard similarity, and gain/loss between repeated isolates.

The important sequencing of logic is this: first we measure the bacterial VF profile, then we compare repeated profiles within residents, and only then do we overlay clinical status such as UTI versus Not_UTI.

**Emphasise:**
- One clinical episode becomes one bacterial genome profile.
- The project is not only status comparison; it is longitudinal VF profiling.
- Clinical status enters as an interpretation layer after the genomic/VF profile is built.

**Avoid:** Do not dive into script names or denominator changes on this slide; keep it plain-English.

**Transition:** Now that the loop is clear, the next slide explains why the numbers change as data move through that loop.

### Slide 3. What changes as data move through the pipeline (1:00)
**Primary message:** Explain that each denominator refers to a different unit of analysis.
**Visual walkthrough:** Use the ladder from clinical episodes to VF-ready episodes, feature space, longitudinal pairs, and clinical transitions.

**Speaker script:**

This slide is here to stop a very common misunderstanding. The counts are not meant to match at every step, because the unit changes.

At the clinical layer we have 583 episodes. Once we require selected assemblies and VF-ready evidence, the analysis dataset has 556 VF-ready episodes. Then we switch from rows to features: 227 VF gene columns are grouped into 32 curated modules.

For the longitudinal analysis, the unit changes again. We are no longer counting episodes; we are counting consecutive within-resident comparisons, which gives 394 longitudinal pairs. Finally, the clinical transition application focuses on a much smaller set: 11 Not_UTI to UTI transition cases.

So the live rule is: always ask what the row is. Is it a clinical episode, an assembly, a VF feature, a pair, or a transition?

**Emphasise:**
- Counts change because the unit changes.
- The same project contains rows, columns, pairs, and transition cases.
- Do not compare denominators as if they are all the same object.

**Avoid:** Do not frame denominator loss as unexplained missingness; much of it reflects unit changes and eligibility filters.

**Transition:** Before we return to VF stability, I want to define what UTI means in this pipeline.

### Slide 4. How an episode becomes UTI in this project (1:00)
**Primary message:** Define the clinical annotation simply: culture support plus catheter-aware compatible symptoms.
**Visual walkthrough:** Show the two criteria feeding into UTI, then the denominator boxes at the bottom.

**Speaker script:**

In this project, UTI is not just bacteria in urine. The primary status requires two things: culture support and compatible symptoms.

Culture support means the urine culture evidence is compatible with possible infection under the primary lower-threshold rule, using >=10^3 CFU/mL support where CFU data are available.

The symptom rule is catheter-aware. Non-catheter and catheter episodes use different clinical logic, because the symptom picture and diagnostic confidence differ between those contexts.

If both culture support and the compatible symptom rule are met, the episode is labelled UTI. Otherwise it is Not_UTI. That gives the primary clinical denominator: 583 clinical episodes, with 18 UTI and 565 Not_UTI.

The key caution is that Not_UTI is a mixed comparator. It can include bacteriuria without the symptom rule, culture-negative or below-threshold episodes, and near-miss or sensitivity-only rows. So Not_UTI should not be described as one uniform biological state.

**Emphasise:**
- UTI = culture support plus compatible clinical picture.
- The symptom rule is catheter-aware.
- Not_UTI is heterogeneous and should be interpreted carefully.

**Avoid:** Do not call every culture-positive episode a UTI, and do not treat Not_UTI as equivalent to healthy/no bacteria.

**Transition:** With the clinical annotation defined, we can return to the main VF-first question: what do repeated urinary E. coli profiles do over time?

### Slide 5. Longitudinal urinary E. coli VF analysis asks about repertoire and stability (1:30)
**Primary message:** Longitudinal VF analysis separates persistent carriage, replacement, and gene-content change.
**Visual walkthrough:** Walk left to right: repeated isolates, three longitudinal patterns, clinical overlay as the later layer.

**Speaker script:**

For each isolate, the VF profile is a measured gene-content snapshot. If we only compare UTI and Not_UTI episodes cross-sectionally, we miss the temporal biology: a resident may carry a stable E. coli profile, may acquire a replacement lineage, or may show apparent VF gene gains or losses between visits.

The middle column names the three patterns I want the listener to keep in mind: persistent lineage, replacement, and gene-content change. These are interpretive buckets, not mutually perfect mechanistic categories.

The right-hand panel is deliberately labelled clinical overlay comes later. UTI_Status is useful, but in this analysis it is a cautious annotation layer placed on top of the genomic/VF pipeline rather than the whole thesis.

**Emphasise:**
- VF presence/absence is a genomic feature, not expression or activity.
- The longitudinal structure is the reason the dataset is interesting.

**Avoid:** Do not say stable VF equals asymptomatic carriage or replacement equals new infection without supporting genomic context.

**Transition:** The next slide shows what data asset feeds that longitudinal question.

### Slide 6. VF-ready data link repeated urinary E. coli genomes to clinical episodes (1:30)
**Primary message:** Define the unit of analysis and the VF-ready denominator.
**Visual walkthrough:** Use the metric tiles first, then the repeated isolate timeline.

**Speaker script:**

The analysis dataset is episode-level. Each row represents a selected urinary E. coli assembly linked to a clinical episode and translated into a binary VF profile.

The VF-ready dataset contains 556 episode-level selected assemblies from 162 participants. Those assemblies are represented by 227 binary VF gene columns. The curation framework groups those genes into 32 modules, including 18 UPEC-candidate modules.

The timeline underneath is schematic: one resident can contribute repeated Not_UTI and UTI-labelled episodes. The key point is that status labels can be overlaid on repeated isolates, but the molecular unit remains the episode-level isolate profile.

**Emphasise:**
- The 162 count is the VF-ready dataset participant count.
- Avoid displaying or discussing the unresolved total clinical participant count.

**Avoid:** Do not imply that all 583 clinical episodes have VF-ready WGS evidence.

**Transition:** Now I’ll go through the exact handoffs that take us from clinical classification to VF outputs.

### Slide 7. Numbered clinical-to-VF pipeline with visible denominators (2:30)
**Primary message:** Orient the colleague to each pipeline handoff and why counts change.
**Visual walkthrough:** Read the flow in numbered order; pause on the places where the denominator changes.

**Speaker script:**

This is the slide I would expect a new colleague to come back to. It shows the current pipeline with the count at each handoff, so they can understand what is being compared.

We begin with 585 classified clinical episodes before primary manual exclusions. After primary clinical inclusion, 583 episodes remain: 18 UTI and 565 Not_UTI, with two excluded. At the WGS QC stage there are 1,291 assembly-level QC records, which is larger because assembler alternatives are still represented.

Canonical selection brings this to 556 selected episode-level assemblies. Those become 556 VF presence/absence rows, and then 556 VF/model-ready episodes. In that model-ready clinical overlay there are 17 UTI and 539 Not_UTI episodes, meaning 27 clinical episodes do not have VF-ready evidence.

The feature framework maps 227 VF gene columns into 32 modules, including 18 UPEC-candidate modules. The longitudinal layer then uses 394 consecutive within-resident comparisons from 144 participants. Finally, the clinical transition application identifies 11 Not_UTI to UTI clinical transitions, 10 of which are WGS/VF-linked, with one missing endpoint.

**Emphasise:**
- Different rows mean different units: clinical episodes, assembly QC records, selected assemblies, VF rows, longitudinal pairs.
- The main deck should not show disputed participant totals or ST totals.

**Avoid:** Do not describe 1,291 assembly QC records as 1,291 unique isolates.

**Transition:** Before we interpret VF figures, I want to define what the VF features are and what they are not.

### Slide 8. VF features: genes, curated modules, and longitudinal similarity (2:00)
**Primary message:** Separate raw gene detection, curated modules, and downstream outputs.
**Visual walkthrough:** Follow the three main boxes from gene matrix to modules to outputs, then explain the caution box.

**Speaker script:**

The raw VF representation is a binary presence/absence matrix: for each selected assembly, each of 227 VF gene columns is marked detected or not detected. That is a genomic detection framework, not a direct measure of transcription, protein production, or activity.

The curated module layer helps make the feature space biologically interpretable. Rather than discussing hundreds of genes one by one, we can summarise systems such as adhesion, iron acquisition, secretion, toxin, capsule or other surface-related features, and unassigned genes.

The outputs are burden, prevalence, Jaccard similarity, gain/loss, and exploratory models. The caution is important: 25.1 percent of the VF matrix is unassigned, so total burden should be interpreted separately from curated module or UPEC-candidate summaries.

**Emphasise:**
- Modules are curation units, not validated disease-causality scores.
- Jaccard similarity measures profile overlap, not strain identity by itself.

**Avoid:** Do not let the audience hear 'module' as 'clinically validated pathway'.

**Transition:** With the feature framework defined, the first evidence slide simply asks what VF genes are common.

### Slide 9. Common VF genes describe the observed repertoire, not disease association (1:30)
**Primary message:** The top-gene plot is descriptive repertoire evidence.
**Visual walkthrough:** Point to the ranked bars and the descriptive label.

**Speaker script:**

This figure ranks the most prevalent VF genes among the 556 VF/WGS-linked urinary E. coli isolates. It is a repertoire view: what genes are commonly detected in this dataset?

You can use it to orient the audience to the feature space: common fimbrial or adhesion-associated genes, iron acquisition systems, and other VFDB-derived features appear across many isolates. But this is not an association result. A gene being common does not mean it distinguishes UTI from Not_UTI or causes symptoms.

If someone asks why this matters, the answer is that it defines the background VF landscape before we start discussing longitudinal stability or clinical overlays.

**Emphasise:**
- Descriptive prevalence comes before association testing.
- The y-axis ranking is about frequency in the VF-ready dataset.

**Avoid:** Do not call any prevalent gene a UTI marker from this slide.

**Transition:** The next slide groups those individual genes into biological modules to make interpretation easier.

### Slide 10. Curated modules organise genes into interpretable biological systems (1:30)
**Primary message:** Modules reduce complexity while preserving biological interpretability.
**Visual walkthrough:** Use the bar heights to show which module categories contain more genes.

**Speaker script:**

This slide shows the module framework. The goal is not to replace gene-level information, but to give the analysis a more interpretable vocabulary.

A new colleague can use this figure to understand which biological systems are represented in the VF matrix. Some modules contain many genes, so they will naturally contribute more to burden summaries, while smaller modules may be biologically meaningful but numerically less dominant.

The important phrasing is curated framework. These modules help us organise the analysis; they are not validated predictors of symptomatic infection.

**Emphasise:**
- Modules are helpful for biological discussion and handover.
- Gene count per module affects how burden-like summaries behave.

**Avoid:** Do not compare module size as though it directly measures importance.

**Transition:** Once the repertoire and modules are defined, the central longitudinal question is how much they change over time.

### Slide 11. Most repeated within-resident VF profiles are highly stable (2:00)
**Primary message:** The central longitudinal result: VF profiles are often conserved within residents.
**Visual walkthrough:** Point to the mass of comparisons near Jaccard 1.000 and the central finding box.

**Speaker script:**

This is one of the most important slides in the deck. Across 394 consecutive within-resident comparisons from 144 participants, the median within-resident VF Jaccard similarity is 1.000. In other words, for the typical consecutive pair, the detected VF gene set is identical.

The additional headline is that 62.4 percent of consecutive comparisons show no VF change. That tells us that the dominant pattern in this dataset is profile stability, not repeated major VF gain or loss.

The careful interpretation is that VF stability supports persistent-profile explanations, but it does not prove same-strain persistence on its own. For that, we need SNP distance, sequence type, and broader genomic context.

**Emphasise:**
- Median Jaccard similarity: 1.000.
- 62.4 percent with no VF change.
- This is descriptive longitudinal evidence.

**Avoid:** Do not claim that stable VF profile proves the same strain persisted.

**Transition:** The next slide asks where change does appear and how we should treat it.

### Slide 12. When profiles change, gain/loss summaries flag candidates for follow-up (2:00)
**Primary message:** VF change events are leads, not automatic biological conclusions.
**Visual walkthrough:** Explain the gain/loss distributions and point to outlier behaviour.

**Speaker script:**

This figure focuses on consecutive pairs where the VF profile does change. It summarises apparent gene gains and losses between repeated isolates.

Most comparisons are close to zero, consistent with the stability shown on the previous slide. But the tails matter. Pairs with many apparent gains or losses are candidates for follow-up, because they may represent strain replacement, changes in assembly or gene-calling confidence, or true gene-content change.

The correct language here is flagging and prioritisation. This slide helps us identify which cases deserve closer review; it does not by itself adjudicate mechanism.

**Emphasise:**
- Gain/loss summaries should trigger context review.
- Outliers are useful even if they are not immediately mechanistic.

**Avoid:** Do not call every gain/loss event within-host evolution.

**Transition:** To interpret those changes, we need lineage context.

### Slide 13. Sequence-type consistency helps interpret VF stability (1:30)
**Primary message:** Lineage context helps distinguish persistence-like from replacement-like comparisons.
**Visual walkthrough:** Contrast same-ST and different-ST comparisons without overreading the plot.

**Speaker script:**

This diagnostic plot asks whether VF similarity behaves differently when consecutive isolates share the same sequence type versus when they do not. It provides lineage context for the stability result.

The broad expectation is that same-ST pairs should often have more similar VF profiles, while different-ST pairs are more likely to reflect replacement-like events or broader gene-content differences. But ST is still a coarse label. It helps triage interpretation; it does not prove persistence alone.

This is why the pipeline combines profile similarity with SNP distance and case-level interpretation when discussing transitions.

**Emphasise:**
- ST is a context variable, not a final answer.
- This slide prepares the audience for cautious clinical interpretation.

**Avoid:** Do not equate same ST with same strain.

**Transition:** Only after the VF and lineage context are established do I introduce clinical-status comparison.

### Slide 14. UTI status overlay remains exploratory under the sparse denominator (2:00)
**Primary message:** No global VF association remains significant after FDR correction.
**Visual walkthrough:** Explain nominal screening versus participant-aware model evidence; point to the safe-claim box.

**Speaker script:**

This is the first clinical-status application. The model-ready clinical overlay contains 556 episodes: 17 UTI and 539 Not_UTI. That imbalance is the central reason for caution.

The figure compares exploratory gene screening with participant-aware modelling. Some genes may look interesting nominally, but the current safe conclusion is that no global VF association remains significant after FDR correction.

So this slide should be used for hypothesis generation and for explaining uncertainty. It is not a claim that the pipeline has found a robust VF marker of symptomatic UTI.

**Emphasise:**
- 17 UTI versus 539 Not_UTI in the VF/model-ready dataset.
- No global VF association survives FDR correction.
- Clinical overlay is exploratory.

**Avoid:** Do not mention older ASB-versus-UTI conclusions or significant lpf claims as current results.

**Transition:** The next slide makes the clinical overlay concrete with one participant transition.

### Slide 15. Participant 20026 illustrates stable VF profile despite symptom emergence (2:00)
**Primary message:** A concrete transition shows symptoms can emerge without detectable VF gain.
**Visual walkthrough:** Walk across the timeline, then the four green evidence tiles, then the two summary boxes.

**Speaker script:**

Participant 20026 is a useful worked example because it links the abstract pipeline to a real longitudinal transition. The participant moves from Not_UTI to UTI over 42 days.

The genomic and VF evidence shows 5 SNPs and a stable VF profile: zero VF genes gained and zero lost, with stable VF modules and strong same-strain evidence. That pattern supports the idea that symptom emergence does not always require a detectable VF repertoire change.

At the casebook level, there are 11 clinical Not_UTI to UTI transitions, 10 WGS/VF-linked transitions, and one missing endpoint. Four transitions show same-strain stable profiles, while three are consistent with strain replacement. These buckets organise evidence. They do not prove mechanism.

The careful interpretation is that low SNP distance plus stable VF profile supports host-state, expression, regulation, inoculum, or other unmeasured explanations, but it does not prove any one of them.

**Emphasise:**
- Not_UTI to UTI over 42 days.
- 5 SNPs, stable VF profile, zero gains/losses.
- Support, not proof, of host-state or unmeasured-regulation hypotheses.

**Avoid:** Do not say stable VF proves host mechanism or symptom causality.

**Transition:** I’ll close the spoken section by separating what is established from what remains uncertain.

### Slide 16. What this VF-first review establishes and what remains uncertain (1:00)
**Primary message:** Close with evidence, boundary, and next steps.
**Visual walkthrough:** Use the three columns as the closing structure.

**Speaker script:**

The evidence column is the firm ground: the current pipeline produces a coherent VF-ready longitudinal dataset with 556 isolate profiles, 227 VF gene columns, 32 modules, and 394 consecutive within-resident comparisons.

The boundary column is equally important. VF presence/absence does not measure expression or activity, and clinical-status association remains exploratory because the UTI denominator is sparse.

The next-step column is where a new colleague can contribute: lineage-aware longitudinal follow-up, expression or regulation hypotheses, and targeted review of transition case studies. The main discussion question is what additional evidence would distinguish stable carriage with host-state change from unmeasured bacterial regulation or replacement.

**Emphasise:**
- End on usable pipeline plus honest uncertainty.
- Invite discussion around next evidence, not around overclaiming current associations.

**Avoid:** Do not close as though the VF analysis has identified definitive UTI drivers.

**Transition:** From here, the appendix is available for operational handover and detailed questions.

## Appendix / Q&A Slides

### Slide 17. Practical map: where to enter and rerun the VF pipeline
**When to open:** Open when the colleague asks how to navigate the project or reproduce outputs.
**Safe claim:** These are current navigation anchors and rerun order, not new analysis results.
**Avoid:** Do not ask the colleague to start from old ASB scripts or legacy figures.

**Speaker script:**

This is the handover slide. I would use it after the main talk, not during the core narrative unless the audience asks about reproducibility.

The key message is that the colleague should start from the canonical VF dataset and diagnostics, then use the pipeline documentation and final figure pack to trace outputs. The final figure pack script is the entry point for current final figures; the VF and longitudinal outputs are the sources of truth for pipeline-specific questions.

### Slide 18. Clinical phenotype denominator and definition
**When to open:** Open when someone asks how UTI and Not_UTI were defined.
**Safe claim:** UTI_Status is current UTI versus Not_UTI; it is not the older ASB-versus-UTI storyline.
**Avoid:** Do not treat Not_UTI as a single biologically uniform state.

**Speaker script:**

This appendix slide explains the clinical annotation layer. It shows the validated clinical denominator: 583 primary clinical episodes, with 18 UTI and 565 Not_UTI.

The clinical rule is catheter-aware and symptom-supported, with culture support. The important framing is that Not_UTI is heterogeneous: it can include culture-supported bacteriuria without compatible symptoms and other non-UTI states. In this VF-first deck, this status layer is used later as an exploratory annotation.

### Slide 19. Module prevalence by status
**When to open:** Open for questions about whether module prevalence differs by UTI_Status.
**Safe claim:** Exploratory clinical annotation only.
**Avoid:** Do not call module differences validated UTI signatures.

**Speaker script:**

This plot gives module prevalence by clinical annotation. It is useful for visual pattern recognition, but the label is exploratory for a reason.

The safe interpretation is that modules can be compared descriptively across UTI and Not_UTI episodes, but sparse UTI counts, repeated measures, and lineage structure prevent strong clinical conclusions from this slide alone.

### Slide 20. Full transition mechanism casebook
**When to open:** Open when someone wants to inspect all Not_UTI to UTI transition categories.
**Safe claim:** 11 clinical transitions, 10 WGS/VF-linked, with stable-profile and replacement-consistent categories.
**Avoid:** Do not imply that every transition has a single settled biological mechanism.

**Speaker script:**

This is the aggregate transition casebook. It summarises clinical Not_UTI to UTI transitions and the genomic/VF evidence attached to them.

Use it to show that transitions are heterogeneous. Some are consistent with same-strain stable-profile transitions, some with replacement, and some remain less certain. The purpose is evidence organisation, not mechanism proof.

### Slide 21. Stable strain and changing clinical state
**When to open:** Open when someone asks how low SNP distance and stable VF relate to symptom emergence.
**Safe claim:** Stable measured profiles support alternative hypotheses; they do not prove them.
**Avoid:** Do not say host state is proven.

**Speaker script:**

This slide expands the idea shown in participant 20026. It illustrates cases where low SNP distance and stable measured VF profile coexist with changing clinical state.

The interpretive value is that symptoms can change even when the measured VF repertoire does not. That points to host context, regulation, expression, inoculum, tissue state, or other unmeasured factors as plausible follow-up hypotheses.

### Slide 22. Population-level robustness boundary
**When to open:** Open for statistical robustness questions about the global VF overlay.
**Safe claim:** No adjusted VF association is confirmatory with only 17 VF-ready UTI episodes.
**Avoid:** Do not overinterpret nominal effects.

**Speaker script:**

This slide reinforces the population-level boundary. The analysis is underpowered for robust clinical-status association because there are only 17 VF-ready UTI episodes.

Use the slide to explain why the deck presents population-level VF findings as exploratory and why no global VF signal should be treated as confirmatory after FDR correction.

### Slide 23. Lineage structure is an interpretation check
**When to open:** Open when someone asks about population structure or confounding by lineage.
**Safe claim:** Diagnostic/descriptive lineage context.
**Avoid:** Do not use the PCoA as a disease-association result.

**Speaker script:**

This PCoA is a diagnostic view of VF profile structure by sequence type. It reminds us that VF profiles are not randomly distributed across isolates; they are shaped by lineage.

That matters because a clinical-status comparison can be confounded by lineage composition. The correct use is diagnostic: lineage structure needs to be considered before clinical interpretation.

### Slide 24. AMR backup only: exploratory transition-level ResFinder context
**When to open:** Open only if asked about AMR or accessory/plasmid changes.
**Safe claim:** Exploratory transition-level context only.
**Avoid:** Do not claim global AMR association or integrate AMR as a primary spoken result.

**Speaker script:**

This is intentionally an appendix-only backup. ResFinder and accessory/plasmid context can help describe individual transitions, but the current deck does not present a completed global VF-plus-AMR association analysis.

If asked, say that AMR is useful context for transition-level interpretation, but it should not be elevated to a primary result in this VF pipeline review.

### Slide 25. Sources of truth and lookup resources for questions
**When to open:** Open at the end or during Q&A when someone asks where a number, figure, or caveat comes from.
**Safe claim:** Use these files to verify provenance and avoid stale claims.
**Avoid:** Do not use backup diagnostics to create new claims during the presentation.

**Speaker script:**

This is the evidence registry. It points the colleague to the count validation files, denominator flow, VF diagnostics, WGS registry, longitudinal summaries, module notes, casebook, and phenotype explanation.

The most important handover message is that questions should be answered from current sources of truth rather than from older exploratory slides or superseded analyses.
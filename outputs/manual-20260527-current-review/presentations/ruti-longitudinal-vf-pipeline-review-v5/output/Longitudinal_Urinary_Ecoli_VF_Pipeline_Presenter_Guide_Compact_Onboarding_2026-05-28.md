# Longitudinal urinary E. coli VF pipeline review: compact onboarding presenter guide

**Purpose:** Help introduce a new team member to the current VF-first longitudinal pipeline and its cautious clinical annotation layer.

**Timing:** 20-22 min presentation plus discussion.

## Avoid Saying
- Do not frame the current story as an ASB-versus-UTI analysis.
- Do not claim a significant lpf result as current.
- Do not claim a global AMR association.
- Do not say stable VF proves host mechanism.

## Slide 1. Virulence-factor profiling of longitudinal urinary E. coli isolates
**Target time:** 1 min

**Primary message:** Open with a VF-first project frame: the aim is to understand the bacterial VF repertoire and its longitudinal stability, with clinical status used later as an annotation.

**How to say it live:**
I would start by saying: this project is not only a recurrent UTI story. The core object is a longitudinal set of urinary E. coli isolates, each converted into a virulence-factor profile. The three questions are: what VF repertoire is present, how stable those profiles are across repeated isolates from the same resident, and how carefully we can overlay UTI versus Not_UTI status. The status comparison matters, but it enters after we understand the data product and the longitudinal bacterial signal.

**Transition:** Before looking at plots, I want to give the plain-English loop of what the pipeline is doing.

**Caution:** Do not let the title drift back into a recurrent-UTI-only framing.

## Slide 2. The YELLOW routine in one plain-English loop
**Target time:** 1.5 min

**Primary message:** Give the newcomer a mental model for the entire routine without scripts or dense counts.

**How to say it live:**
Talk through the slide left to right. A clinical episode is a urine-sampling event. From that episode, we have a urinary E. coli isolate. Sequencing turns that isolate into a genome assembly. The VF pipeline asks which virulence-factor genes are detected in that assembly, so one isolate becomes one VF gene row. Those rows can then be summarised into modules, burden, similarity, and gain or loss. Only after that do we overlay clinical status. The key phrase is: one episode becomes one bacterial profile, and repeated profiles are compared over time.

**Transition:** The next slide explains why the denominator changes depending on which object we are counting.

**Caution:** Do not mention script names here; this slide is intentionally non-technical.

## Slide 3. Unit changes explain denominator changes
**Target time:** 1.5 min

**Primary message:** Make denominator changes feel logical rather than alarming.

**How to say it live:**
This slide is useful for preventing confusion. The project begins with 583 included clinical episodes, but not every clinical episode has a selected VF-ready assembly. That is why the VF-ready dataset is 556 episodes. From those 556 rows, we represent 227 VF genes and 32 modules. When we move to longitudinal analysis, the unit changes again: now we are counting consecutive within-resident pairs, giving 394 comparisons. When we focus on Not_UTI to UTI clinical transitions, the unit narrows to 11 transition cases. So the counts are not contradictory; they answer different questions.

**Transition:** Because UTI status is used later, I will define exactly how an episode receives that label.

**Caution:** Avoid comparing the 583, 556, 394, and 11 denominators as if they are the same denominator.

## Slide 4. How an episode becomes UTI in this project
**Target time:** 1.5 min

**Primary message:** Clarify that UTI is defined by culture support plus catheter-aware symptoms.

**How to say it live:**
Here I would say: in this analysis, UTI is not simply bacteria in urine. A UTI episode requires culture support and compatible symptoms, and the symptom rule is catheter-aware. If both conditions are met, the episode is UTI. If not, it sits in Not_UTI. That means Not_UTI is deliberately heterogeneous: it can include bacteriuria without the symptom rule, episodes that do not meet culture support, or other near-miss contexts. This matters because later status comparisons are comparing a small UTI group against a broad comparator.

**Transition:** Now that the clinical annotation is defined, we can return to the VF-ready data asset.

**Caution:** Do not describe Not_UTI as healthy controls.

## Slide 5. VF-ready data asset
**Target time:** 1 min

**Primary message:** Show what the dataset physically represents: repeated isolates with selected assemblies and VF profiles.

**How to say it live:**
The working VF dataset has 556 episode-level selected assemblies from 162 participants. Each row is a selected assembly linked back to the clinical episode. Across those rows, the pipeline records 227 binary VF gene columns, then groups those genes into 32 curated modules, including 18 UPEC-candidate modules. The timeline schematic is deliberately simple: repeated isolates from the same resident can carry different status labels, but the bacterial profiles are what we compare longitudinally.

**Transition:** The next slide explains how those VF profiles are encoded and summarised.

**Caution:** Do not display or state a total clinical participant count outside the VF-ready dataset.

## Slide 6. VF features: genes, modules, and similarity
**Target time:** 1.5 min

**Primary message:** Explain the feature representation and why modules help interpretation.

**How to say it live:**
This is the translation from genome to analysis features. The simplest object is the binary gene matrix: present or absent for each VF gene in each isolate. Because individual genes can be hard to interpret, the pipeline also maps genes into curated biological modules. Those modules help us talk about adhesion, iron acquisition, secretion, capsule and surface structures, toxins, and unassigned material. Then the longitudinal outputs compare profiles across time using burden, prevalence, Jaccard similarity, and gain/loss summaries.

**Transition:** With that feature language in place, we can look at the observed VF repertoire.

**Caution:** Modules are curation units, not validated disease prediction scores.

## Slide 7. VF repertoire summarised as genes and modules
**Target time:** 1.5 min

**Primary message:** Separate descriptive prevalence from association testing.

**How to say it live:**
Use the left panel as a descriptive prevalence ranking: these are common VF genes among the 556 VF-ready isolates. That does not mean they cause symptoms or distinguish UTI. The right panel shows how genes are distributed across curated modules. I would emphasise that this gives the colleague vocabulary for later slides. We are first learning what is present and how it is organised, before asking whether any clinical annotation appears to line up with it.

**Transition:** The central result comes next: how similar repeated profiles are within the same resident.

**Caution:** Do not interpret high prevalence as clinical importance by itself.

## Slide 8. Most repeated within-resident VF profiles are highly stable
**Target time:** 2 min

**Primary message:** Present the main longitudinal VF result.

**How to say it live:**
This is one of the main biological messages. The analysis considers 394 consecutive within-resident comparisons from 144 participants. The median Jaccard similarity is 1.000, and 62.4% of consecutive comparisons show no VF change. In plain terms, many repeated urinary E. coli isolates from the same resident have extremely similar VF presence/absence profiles. That suggests that, for many residents, the measured VF repertoire is conserved over time.

**Transition:** The next question is what we do with the comparisons that are not perfectly stable.

**Caution:** High VF similarity supports stability, but does not prove same-strain persistence without SNP or lineage context.

## Slide 9. Gain/loss summaries flag candidates for follow-up
**Target time:** 1.5 min

**Primary message:** Explain how observed VF changes should be interpreted cautiously.

**How to say it live:**
This plot looks at where VF genes appear to be gained or lost across consecutive pairs. I would frame this as a triage view, not a final mechanism. A gain or loss can happen because the resident has a replacement isolate, because of assembly or calling variation, or because there is genuine gene-content change. The value of the slide is that it tells us where to look more closely, especially when combined with ST and SNP context.

**Transition:** That is why we next bring in sequence-type consistency.

**Caution:** Avoid saying that a plotted gain or loss automatically represents biological evolution in the resident.

## Slide 10. Sequence-type consistency helps interpret VF stability
**Target time:** 1.5 min

**Primary message:** Use lineage as context for interpreting profile stability or change.

**How to say it live:**
Here the point is that VF profile stability needs lineage context. If repeated isolates share the same sequence type and have very similar VF profiles, that supports a persistent-lineage interpretation. If sequence type changes, then a VF change may be more consistent with replacement. But ST is still a coarse diagnostic. It helps organise the evidence; it does not independently prove same strain.

**Transition:** Only after this VF and lineage structure do we move to the clinical-status overlay.

**Caution:** Do not describe ST agreement as definitive persistence.

## Slide 11. Clinical-status VF signals remain exploratory
**Target time:** 2 min

**Primary message:** Make the status analysis useful but bounded.

**How to say it live:**
This is where UTI versus Not_UTI enters the analysis as a clinical annotation. The VF/model-ready dataset has 556 episodes, but only 17 are UTI and 539 are Not_UTI. The three blocks show the intended logic more clearly than the dense model plot: first, a simple screen can produce nominal-looking VF differences; second, after participant-aware modelling and FDR correction, there are 0 robust global VF hits; third, sparse/separation flags and lineage structure mean we should treat this as a prioritisation layer. The safe conclusion is that the current data support hypothesis generation, not a confirmed UTI-associated VF signature.

**Transition:** The next slide shows why the longitudinal case view is still biologically useful.

**Caution:** Do not mention old significant lpf claims or imply a corrected global association.

## Slide 12. Participant 20026 transition example
**Target time:** 2 min

**Primary message:** Make the longitudinal argument concrete.

**How to say it live:**
This is a worked example. Participant 20026 moves from Not_UTI to UTI over 42 days. The isolates are separated by only 5 SNPs and show no VF gene gain or loss in the measured profile. That means the measured virulence repertoire appears stable while the clinical state changes. The careful interpretation is that stable VF plus low SNP distance supports a host-state, timing, or unmeasured-regulation hypothesis. It does not prove the mechanism. Across the broader casebook there are 11 clinical Not_UTI to UTI transitions, 10 linked to WGS/VF evidence, with 4 stable-profile and 3 replacement-consistent transitions.

**Transition:** I would close the spoken deck by turning those observations into takeaways and open decisions.

**Caution:** Never say stable VF proves host mechanism.

## Slide 13. Takeaways, limitations, and next steps
**Target time:** 1 min

**Primary message:** Close the spoken narrative with three clear conclusions.

**How to say it live:**
The wrap-up is: first, the current pipeline produces a usable longitudinal VF dataset. Second, VF profiles are commonly stable within residents across repeated urinary E. coli isolates. Third, clinical-status associations remain exploratory because the UTI denominator is sparse and repeated measures and lineage structure matter. The next research step is not to overclaim a single VF marker, but to use this pipeline to prioritise lineage-aware follow-up, expression or regulation hypotheses, and careful review of transition cases.

**Transition:** Then pause for questions, using the appendix slides depending on what the colleague asks.

**Caution:** Do not turn the limitations into an apology; they define the correct next analyses.

## Slide 14. Appendix: practical handover map
**Target time:** Q&A

**Primary message:** Use when the colleague asks where to start or how to rerun the project.

**How to say it live:**
Point them first to the canonical outputs: status_map.csv for clinical annotation, vf_pa_all.csv for the raw VF matrix, vf_analysis_ready.csv for the 556-row analysis dataset, and the gene/module map for curation. For reruns, walk from clinical status through assembly selection, VF matrix building, analysis-ready construction, modules and scores, longitudinal summaries, then the final figure pack.

**Transition:** Use this slide as the practical handover rather than trying to explain the repo from memory.

**Caution:** Do not ask a new colleague to begin with every script; begin with outputs and then trace backward.

## Slide 15. Appendix: detailed clinical-to-VF pipeline
**Target time:** Q&A

**Primary message:** Use for behind-the-scenes denominator and script-order questions.

**How to say it live:**
This is the dense map. It shows why the simplified count ladder is true: 585 classified episodes before exclusions, 583 primary included clinical episodes, 1,291 assembly-level QC records including assembler alternatives, 556 selected assemblies and VF rows, then the 556-row VF/model-ready dataset. It also shows the feature framework and the transition application. This is the slide to use when someone asks exactly where rows are lost or why assembly records exceed episode records.

**Transition:** If they ask about transition mechanisms specifically, move to the casebook slide.

**Caution:** Keep the unit-of-analysis language explicit.

## Slide 16. Appendix: full transition mechanism casebook
**Target time:** Q&A

**Primary message:** Use for detailed questions about Not_UTI to UTI transitions.

**How to say it live:**
This appendix slide shows the broader transition casebook. The safe statement is that 11 clinical Not_UTI to UTI transitions were identified; 10 have WGS/VF-linked evidence and 1 is missing the endpoint. Four transitions are stable-profile cases and three are replacement-consistent. These are evidence buckets, not proven mechanisms.

**Transition:** For population-level signal or confounding questions, go to the final appendix slide.

**Caution:** Do not make the transition buckets sound mutually exhaustive or mechanistically proven.

## Slide 17. Appendix: robustness and lineage diagnostics
**Target time:** Q&A

**Primary message:** Use when someone asks whether the status signal is robust or lineage-confounded.

**How to say it live:**
The left panel reinforces the population-level boundary: no global VF association is significant after FDR correction. The right panel is a lineage diagnostic, showing why sequence type and bacterial population structure must be considered before interpreting status differences. This slide is there to keep the discussion honest if someone asks whether there is a global VF signature of UTI in the current dataset.

**Transition:** Return to the main takeaway: strong descriptive VF pipeline, cautious clinical overlay.

**Caution:** Do not introduce a global VF-plus-AMR claim; AMR is not part of this compact spoken deck.

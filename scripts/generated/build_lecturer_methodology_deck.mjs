import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

process.on("uncaughtException", (error) => {
  console.error(`BUILD_ERROR: ${error?.message ?? error}`);
  process.exit(1);
});
process.on("unhandledRejection", (error) => {
  console.error(`BUILD_ERROR: ${error?.message ?? error}`);
  process.exit(1);
});

const tmp = "/var/folders/fp/wwwk1rbj70l0k92s92kd2z500000gp/T/codex-presentations/manual-20260710-lecturer-methods/tmp";
const starter = path.join(tmp, "template-starter.pptx");
const layoutDir = path.join(tmp, "template-starter-layout");
const previewDir = path.join(tmp, "final-preview");
const finalLayoutDir = path.join(tmp, "layout", "final");
const finalPptx = "/Users/Aamir/Desktop/rUTIs/outputs/lecturer_methodology_pack/rUTI_complete_methodology_for_lecturer.pptx";

const content = {
  1: [
    "AUDITED METHODOLOGY / LECTURER DISCUSSION",
    "Audited methods for longitudinal urinary E. coli",
    "One clinical definition, one selected genome per episode, and a separate Longcycler-only sensitivity analysis",
    "Clinical definition",
    "583 episodes classified before genomics",
    "UTI requires culture support plus compatible catheter-aware symptoms.",
    "Executed analysis",
    "556 mixed-canonical genomes",
    "The selected set contains 532 Longcycler assemblies and 24 Flye fallbacks.",
    "Sensitivity analysis",
    "532 Longcycler-selected genomes",
    "This route is rebuilt separately and compared with—not substituted for—the executed mixed analysis.",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "Internal scientific audit | regenerated 10 July 2026 | interpretation remains non-causal",
  ],
  2: [
    "PROJECT ORIENTATION",
    "Different numbers count different things",
    "Before interpreting a result, name the unit, filter and denominator.",
    "Sources: status map, canonical assembly selection, VF matrix and fresh pairwise outputs",
    "02",
    "Research question and unit rule",
    "Do closely related urinary E. coli genomes persist across a resident's sampled episodes? Every answer must name the unit, retained rows, filter and denominator.",
    "Clinical episodes",
    "583",
    "included visits or events with a primary UTI label",
    "Selected genome rows",
    "556",
    "one QC-passing canonical assembly per retained genomic episode",
    "VFDB feature space",
    "227",
    "binary detected gene features measured across selected genomes",
    "Adjacent transitions",
    "394",
    "mixed-pipeline comparisons between adjacent retained episodes",
    "Genome pairs",
    "963",
    "all unordered within-participant canonical genome comparisons",
    "One participant can contribute several episodes, genome pairs and transitions.",
  ],
  3: [
    "CLINICAL CLASSIFICATION",
    "UTI requires culture support plus compatible symptoms",
    "The clinical label is constructed independently of assembler choice.",
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0Source: results/clinical/status_map.csv | primary catheter-adjusted definition",
    "03",
    "1. Culture support",
    "Recorded CFU is used where available; a documented Beoord-category lower bound is used when CFU is absent.",
    ">=10^3 CFU/mL lower-bound rule",
    "+",
    "2. Symptom rule",
    "Non-catheter: local symptoms, or flank pain plus systemic symptoms. Catheter: systemic symptoms.",
    "catheter-aware symptom rule",
    "UTI",
    "Both components met",
    "Primary clinical denominator",
    "585 records -> 583 included from 166 participants: 18 UTI and 565 Not_UTI",
    "Not_UTI is not a healthy control",
    "It includes episodes below threshold or with indeterminate evidence; it must not be treated as a healthy comparator.",
    "Clinical label fixed",
    "Genome linked later",
  ],
  4: [
    "ASSEMBLY SELECTION / DETAILED FLOW",
    "Two evidence tracks converge into one episode-level genomic dataset",
    "Clinical status is assigned independently; assembly alternatives are QC-screened before one genome is selected per participant-timepoint.",
    "Sources: metadata_fasta_discovery_manifest.csv; assembly_metadata.csv; canonical_assembly_selection.csv; status_map.csv",
    "04",
    "\u00A0\u00A0\u00A0\u00A0FASTA inventory",
    "1,303 files discovered, including files later excluded or unmatched",
    "1",
    "00_make_assembly_metadata.r",
    "\u00A0\u00A0\u00A0\u00A0Metadata reconciliation",
    "1,299 metadata rows: 1,295 matched FASTAs plus 4 expected/no-FASTA audit rows",
    "2",
    "assembly_metadata.csv",
    "\u00A0\u00A0\u00A0\u00A0Primary candidates",
    "1,291 candidate assembly records across 579 episode keys",
    "3",
    "manual curation + linkage",
    "\u00A0\u00A0\u00A0\u00A0Assembly QC",
    "1,211 assembly records pass readability, 4-6 Mb, <=200 contigs and N50 >=20 kb",
    "4",
    "12a_wgs_qc.R",
    "\u00A0\u00A0\u00A0\u00A0Canonical selection",
    "556 selected genomes: one QC-passing assembly per retained episode key",
    "5",
    "canonical_assembly_selection.csv",
    "\u00A0\u00A0\u00A0\u00A0Longcycler component",
    "532 selected genomes use the preferred Longcycler assembly",
    "6",
    "assembler == longcycler",
    "\u00A0\u00A0\u00A0\u00A0Flye component",
    "24 selected genomes are Flye fallbacks",
    "7",
    "assembler == flye",
    "\u00A0\u00A0\u00A0\u00A0Executed mixed set",
    "556 rows / 162 participants / 17 UTI / 539 Not_UTI",
    "8",
    "current generated pipeline",
    "\u00A0\u00A0\u00A0\u00A0Sensitivity set",
    "532 rows / 161 participants / 16 UTI / 516 Not_UTI",
    "9",
    "Longcycler-only restriction",
  ],
  5: [
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0ASSEMBLER SENSITIVITY",
    "The executed analysis is mixed; Longcycler-only tests assembler sensitivity",
    "\u00A0\u00A0\u00A0The restriction reduces one source of technical variation, but it can also change who and what remains represented.",
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0Source: selected canonical assembly manifest and Longcycler sensitivity rebuild",
    "05",
    "556",
    "executed mixed-canonical rows",
    "532",
    "Longcycler-selected rows",
    "24",
    "Flye fallback rows excluded",
    "17 -> 16",
    "UTI rows after restriction",
    "9",
    "new Longcycler adjacencies",
    "Why the timeline is rebuilt",
    "Longcycler",
    "retained T0",
    "Longcycler",
    "retained T1",
    "Flye",
    "removed T2",
    "Longcycler",
    "newly adjacent T3",
    "Longcycler",
    "retained UTI event",
    "Interpretation stance",
    "The restriction removes one participant and one UTI-linked row, and it creates nine new Longcycler-to-Longcycler adjacencies. Rebuild the timeline and treat fallback availability as a possible selection mechanism.",
  ],
  6: [
    "GENOMIC ANALYSES",
    "Selected genomes feed parallel analysis branches",
    "Each branch answers a different question; the branches are integrated only where the analysis explicitly joins them.",
    "\u00A0\u00A0Sources: scripts 02, 06, 11, 12b, 12c and 22",
    "06",
    "VFDB screening",
    "Script 02: ABRicate at >=80% identity and >=80% coverage yields 227 binary detected features.",
    "MLST and integration",
    "Script 06 supplies preferred ST context; script 22 joins clinical status, VF, MLST and selected-assembly metadata.",
    "Genome comparison and context",
    "Script 11: dnadiff + Mash. Separate population branches: Parsnp/snp-dists (12b) and Prokka/Panaroo (12c).",
    "All branches start from the selected canonical set; arrows between methods would incorrectly imply a sequence.",
    "VFDB / 02",
    "MLST / 06",
    "join / 22",
    "dnadiff / 11",
    "Parsnp / 12b",
    "Panaroo / 12c",
    "Interpretation boundary",
    "ST is lineage context. dnadiff SNPs are assembly-to-assembly. Parsnp is core-genome. wgMLST allele distances are separate and were not recalculated here.",
    "No genomic branch on its own establishes why a UTI occurred.",
  ],
  7: [
    "LONGITUDINAL RESULT",
    "Longitudinal comparisons are rebuilt after changing the retained episodes",
    "Fresh SNP counts come only from SHA-256-matched current FASTAs and reports.",
    "Sources: vf_longitudinal_transitions.csv and results/sensitivity/longcycler_only",
    "07",
    "394",
    "mixed transitions / 144 participants",
    "138 / 394",
    "mixed transitions at <=25 dnadiff SNPs",
    "371",
    "Longcycler transitions / 139 participants",
    "140 / 371",
    "Longcycler transitions at <=25",
    "5 / 9",
    "Not_UTI-to-UTI transitions at <=25",
    "How the focused result is obtained",
    "Mixed timeline",
    "394 pairs / 144 participants",
    "Apply LC restriction",
    "532 retained episode rows",
    "Rebuild",
    "+9 newly adjacent LC pairs",
    "LC timeline",
    "371 pairs / 139 participants",
    "UTI transition",
    "5 of 9 meet <=25",
    "Safe interpretation",
    "Fresh result: 138/394 mixed and 140/371 Longcycler transitions meet the operational <=25 rule; 5/9 Longcycler Not_UTI-to-UTI transitions do. These are descriptive comparisons, not proof of uninterrupted carriage or causation.",
  ],
  8: [
    "CONCLUSION / SUPERVISOR DECISIONS",
    "The claim-linked methods support cautious genomic description",
    "The analysis is strongest when the executed route, sensitivity route and interpretation limits remain separate.",
    "Source: audited scripts, refreshed outputs and explicit scope boundaries",
    "08",
    "Supported evidence",
    "Defensible denominators; detected VFDB presence/absence; preferred MLST lineage context; fresh hash-verified dnadiff comparisons; rebuilt mixed and Longcycler timelines.",
    "Interpretation boundary",
    "No antibiotic, demographic, host-factor or wgMLST result was reconstructed. <=25 is operational and needs citation or calibration. Repeated observations limit simple tests.",
    "Lecturer decisions",
    "Approve mixed-primary/Longcycler-sensitivity framing; decide how to calibrate <=25; complete missing laboratory, assembler, database and software provenance.",
    "Discussion prompt: is this framing acceptable for the thesis, and what evidence should be required before making a same-strain or causal claim?",
  ],
  9: [
    "APPENDIX / DENOMINATORS",
    "Key denominators keep units and filters visible",
    "\u00A0\u00A0\u00A0These values are verified from current generated files after the cache-safe rerun.",
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0Sources: status_map, canonical selection, VF-ready and sensitivity count tables",
    "09",
    "Rule for auditing a count",
    "Record the source file, filter, denominator, formula and unit. A change in count is interpretable only after all five are stated.",
    "Clinical records",
    "585 -> 583",
    "two primary exclusions; 166 included participants; 18 UTI / 565 Not_UTI",
    "Assembly records",
    "1,291 -> 1,211",
    "candidate assembler alternatives to implemented assembly-QC pass",
    "Selected genomes",
    "556",
    "one canonical QC-passing assembly per eligible genomic episode; 162 participants; 17 UTI / 539 Not_UTI",
    "Longcycler sensitivity",
    "532",
    "selected Longcycler rows; 161 participants; 16 UTI / 516 Not_UTI",
    "Longitudinal transitions",
    "394 -> 371",
    "mixed adjacent pairs to rebuilt Longcycler adjacent pairs",
    "Different units explain the changes; none of these is a participant-only funnel.",
  ],
  10: [
    "APPENDIX / METHODS TABLE",
    "Report implemented thresholds—and criteria not actually applied",
    "This table describes the claim-linked route; missing versions are not inferred.",
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0Sources: scripts 00b, 02, 06, 11, 12a-12c, 22 and 14",
    "10",
    "Implemented thresholds",
    "Tool or branch",
    "Culture support",
    ">=10^3 CFU/mL lower-bound rule",
    "Genome size",
    "4-6 Mb",
    "Assembly structure",
    "<=200 contigs; N50 >=20 kb",
    "Same-strain support",
    "<=25 dnadiff SNPs; operational, uncalibrated",
    "VFDB screen",
    ">=80% identity and >=80% coverage",
    "Provider MLST QC",
    "PercGoodTargets >=95",
    "Canonical inputs",
    "Longcycler preferred; Flye fallback (12a)",
    "VF + clinical join",
    "ABRicate/VFDB (02) then episode join (22)",
    "Lineage context",
    "Provider SeqSphere; labelled local fallback (06)",
    "Pairwise similarity",
    "MUMmer dnadiff + Mash (11 + helper)",
    "Population context",
    "Parsnp/snp-dists (12b); Prokka/Panaroo (12c)",
    "Exploratory model",
    "Participant random intercept with documented GLM fallback (14)",
    "Criteria not applied",
    "Completeness, contamination and read coverage are not part of the implemented 12a assembly-QC screen.",
  ],
  11: [
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0APPENDIX / DISTANCE METHODS",
    "dnadiff SNPs, core-genome SNPs and wgMLST alleles are separate measures",
    "Never convert one measure into another or apply a threshold across methods without validation.",
    "\u00A0\u00A0Sources: scripts 11 and 12b; wgMLST source data were not available for reanalysis",
    "11",
    "Pairwise dnadiff",
    "Compares two complete assemblies. pairwise_metrics.csv TotalSNPs comes from this branch.",
    "Parsnp core genome",
    "Aligns the shared core across the 556 selected genomes, then snp-dists measures that alignment.",
    "wgMLST alleles",
    "Counts differing called loci under a particular scheme and QC process; not recalculated here.",
    "Input and method define the distance. Similar numerical values do not make the outputs equivalent.",
    "assembly pair",
    "path + SHA-256",
    "core alignment",
    "reference-sensitive",
    "scheme loci",
    "scheme-version-sensitive",
    "Operational rule",
    "The project uses <=25 dnadiff SNPs for strong same-strain support. This is not a universal cutoff and requires a defensible citation or dataset-specific calibration.",
    "Directional corroboration is acceptable; numerical interchangeability is not.",
  ],
  12: [
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0APPENDIX / STATISTICAL CAUTION",
    "The exploratory model cannot turn sparse repeated data into causal evidence",
    "Script 14 is included because its results inform the interpretation boundary.",
    "Source: 14_genotype_phenotype_model.R and current model diagnostics",
    "12",
    "Design",
    "The mixed dataset has 556 rows but only 17 UTI outcomes. A participant random intercept addresses repeated episodes in the mixed model.",
    "Selection and fitting",
    "Features at 5-95% prevalence are screened by Fisher test; p<0.10 or the top 50 enter models with timepoint and batch. GLM is the error fallback.",
    "Interpretation",
    "Fallback models lose clustering; sparse/separated cells and singular fits remain. No feature reaches FDR <0.05, so results are exploratory.",
    "Safe wording: describe associations and uncertainty; do not call a feature a risk factor, mechanism or cause of UTI.",
  ],
  13: [
    "APPENDIX / SCRIPT MAP",
    "Direct claim-producing scripts form the traceability core",
    "Other scripts may be mentioned when they clarify population context, quality checks or interpretation.",
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0Source: claim-linked scripts and their current generated outputs",
    "13",
    "\u00A0\u00A0\u00A0\u00A0Clinical definition",
    "00a cleans records; 00b creates the primary UTI_Status",
    "1",
    "input to all claims",
    "\u00A0\u00A0\u00A0\u00A0Assembly linkage",
    "00_make reconciles prepared FASTAs and episode metadata",
    "2",
    "upstream input",
    "\u00A0\u00A0\u00A0\u00A0QC and selection",
    "12a applies implemented QC and selects one canonical assembly",
    "3",
    "defines 556 genomes",
    "\u00A0\u00A0\u00A0\u00A0VFDB screening",
    "02 creates 227 selected-canonical detected features",
    "4",
    "feature input",
    "\u00A0\u00A0\u00A0\u00A0MLST context",
    "06 integrates preferred provider calls and labelled local fallback",
    "5",
    "lineage input",
    "\u00A0\u00A0\u00A0\u00A0Episode integration",
    "22 joins clinical status, VF, MLST and selected-assembly metadata",
    "6",
    "analysis-ready rows",
    "\u00A0\u00A0\u00A0\u00A0Pairwise comparison",
    "11 + helper compare all 963 canonical within-resident genome pairs",
    "7",
    "fresh provenance",
    "\u00A0\u00A0\u00A0\u00A0Mixed longitudinal",
    "24 orders retained episodes and builds 394 adjacent comparisons",
    "8",
    "executed route",
    "\u00A0\u00A0\u00A0\u00A0Longcycler sensitivity",
    "rebuild_longcycler_sensitivity removes Flye rows and rebuilds adjacency",
    "9",
    "sensitivity route",
  ],
  14: [
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0APPENDIX / EVIDENCE STILL NEEDED",
    "State the assumptions and recover missing provenance before submission",
    "The documented workflow starts from prepared clinical exports, metadata and FASTAs—not raw reads.",
    "\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0Source: methodology audit findings and script-level provenance review",
    "14",
    "Assumptions to state",
    "The clinical rule is valid; one selected assembly represents an episode; missing WGS and Flye fallback are not random; adjacent samples do not prove continuous carriage; <=25 is only an operational support rule.",
    "Evidence still needed",
    "Culture-to-colony selection, ONT chemistry, basecaller/read QC/coverage, assembler and polishing commands, database/scheme releases, software versions, code revision and frozen R environment.",
    "Claims to keep out",
    "Antibiotic effects, demographics, host causation, wgMLST reproduction, mutation-location conclusions and causal VF mechanisms are unsupported here.",
    "Before submission: complete provenance, calibrate or cite <=25, assess missing-WGS/assembler selection, and keep the executed mixed analysis distinct from the Longcycler sensitivity analysis.",
  ],
};

const notes = {
  1: "Simple explanation: the analysis first decides what each clinical episode means, then chooses one genome to represent that episode. The 556-genome mixed set is the executed analysis; the 532-genome Longcycler set is a comparison route.\n\nTechnical detail: the executed set contains 532 selected Longcycler assemblies and 24 selected Flye fallbacks. Longcycler-only is rebuilt after restriction and must not be relabelled as the original primary pipeline.\n\nSafe wording: 'The Longcycler restriction tests sensitivity to assembler mixing; it does not prove that Longcycler is universally superior.'",
  2: "Simple explanation: these numbers count different objects, so they are not one shrinking queue.\n\nTechnical detail: a participant may contribute multiple episodes; an episode may have multiple candidate assemblies but only one selected genome; 227 is a number of feature columns; 963 is every unordered within-participant genome pair; 394 counts adjacent retained-episode transitions.\n\nLikely question: 'Why can 963 pairs come from 556 genomes?' Answer: participants with several genomes contribute every within-person pair, not only neighbouring visits.",
  3: "Simple explanation: a UTI label needs both bacterial-growth evidence and compatible symptoms.\n\nTechnical detail: the culture rule uses a >=10^3 CFU/mL lower bound from recorded CFU, with Beoord-category lower bounds when CFU is absent. Non-indwelling episodes qualify through a local urinary symptom or flank pain plus a systemic symptom; indwelling-catheter episodes require systemic symptoms. Unknown or incomplete evidence is Not_UTI with an indeterminate subgroup. Counts use UTI_Status: 585 records, 583 primary-included episodes, 166 participants, 18 UTI and 565 Not_UTI.\n\nSafe wording: Not_UTI is not a healthy or bacteria-free control.",
  4: "Simple explanation: this is two tracks, not one shrinking funnel. Clinical labels are made independently. In parallel, several assembly files can represent one episode until QC and canonical selection leave one genome. The tracks meet only at the participant-timepoint join.\n\nTechnical detail: 1,303 candidate FASTAs reconcile to 1,299 metadata rows (1,295 matched files plus four expected/no-FASTA audit rows), then 1,291 primary candidate assembly records across 579 episode keys. Implemented QC leaves 1,211 passing assembly records across 556 episode keys. Canonical selection keeps one genome per key: 532 Longcycler plus 24 Flye. Linking to clinical status yields 556 mixed rows; 27 of 583 clinical episodes have no selected WGS.\n\nQC boundary: readability, 4-6 Mb, <=200 contigs and N50 >=20 kb. Completeness, contamination and read coverage are not applied here.\n\nLikely question: 'Why is 1,303 to 1,299 not simple attrition?' Answer: unmatched FASTAs are removed while four audit-only expected rows without FASTAs are present in metadata.",
  5: "Simple explanation: using one assembler is like using one reconstruction method throughout, but removing observations can change the sampled group and which visits become neighbours.\n\nTechnical detail: the restriction removes 24 Flye fallback rows, changes 162 to 161 participants and 17 to 16 UTI rows, and creates nine new Longcycler-to-Longcycler adjacencies. Therefore the Longcycler timeline is rebuilt from 532 retained rows rather than filtered from the mixed transition table.\n\nSafe wording: report this as a sensitivity analysis and acknowledge possible selection by fallback availability.",
  6: "Simple explanation: the same selected genomes are examined through several different lenses. One lens does not manufacture the output of another.\n\nTechnical detail: script 02 creates the selected-canonical ABRicate/VFDB matrix; script 06 provides preferred MLST lineage context; script 22 performs the episode-level clinical/VF/MLST/assembly join. Script 11 and its helper provide hash-verified dnadiff and Mash comparisons. Scripts 12b and 12c provide separate Parsnp core-genome and Prokka/Panaroo pangenome context.\n\nSafe wording: ST is lineage context, VF is detected presence/absence and dnadiff is assembly-to-assembly. None alone explains why UTI occurred.",
  7: "Simple explanation: changing the retained visits changes the before-and-after comparisons, so both timelines must be built independently.\n\nTechnical detail: the mixed cohort has 394 transitions from 144 participants; 138 meet TotalSNPs <=25. The Longcycler rebuild has 371 transitions from 139 participants; 140 meet the rule. Of nine Longcycler Not_UTI-to-UTI transitions, five meet it. The extra strict count in the smaller set is possible because nine new adjacencies are created. Values come from current SHA-256-matched reports.\n\nSafe wording: <=25 is an operational support boundary pending citation or calibration. Adjacent sampled episodes do not prove uninterrupted carriage, and continuity does not prove causation.",
  8: "Simple explanation: the analysis can describe what was detected and how genomes compare, but it cannot establish a cause of UTI.\n\nTechnical detail: defensible claims use the clinical definition, selected canonical inputs, detected VF features, preferred MLST context and fresh dnadiff-supported longitudinal comparisons. Antibiotic, demographic, host-factor and wgMLST claims were not reconstructed.\n\nAsk the lecturer to decide: whether the mixed-primary/Longcycler-sensitivity framing is acceptable, how <=25 should be justified, and which missing laboratory/software records can be recovered.",
  9: "Use this slide when a denominator is challenged. State source, filter, unit and formula before interpreting a change. Clinical records, candidate assemblies, selected genomes, participants and transitions are not interchangeable. The mixed selected row count is 556 from 162 participants with 17 UTI and 539 Not_UTI; the Longcycler sensitivity has 532 rows from 161 participants with 16 UTI and 516 Not_UTI.",
  10: "These are the criteria and tools that directly support presented claims. The implemented assembly screen is readability, 4-6 Mb, <=200 contigs and N50 >=20 kb. Completeness, contamination and read coverage are not implemented by 12a and must not be reported as if they were. The <=25 dnadiff boundary is operational and uncalibrated. Versions, commands and database releases still need to be frozen.",
  11: "Simple explanation: three relatedness methods count different things.\n\nTechnical detail: dnadiff compares an assembly pair and the current cache records both paths and SHA-256 hashes. Parsnp aligns a cohort core and is reference/alignment sensitive. wgMLST counts differing called loci under a named scheme and QC process; its source table was unavailable here.\n\nSafe wording: directional agreement can corroborate relatedness, but the numbers and thresholds are not interchangeable.",
  12: "Simple explanation: many rows do not compensate for having only 17 UTI outcomes, especially when residents appear repeatedly.\n\nTechnical detail: script 14 screens features with 5-95% prevalence, uses Fisher p<0.10 or the top 50 for modelling, and fits feature + timepoint + batch with a participant random intercept. A standard GLM is used after fitting errors, losing clustering. Singular fits, sparse/separated cells and selected-family FDR limit inference; no feature reaches FDR <0.05.\n\nSafe wording: hypothesis-generating association, not a risk factor, mechanism or cause.",
  13: "This is claim traceability, not a list of every R file. Read each card as 'claim or dataset -> script that creates it'. The core chain is 00a/00b, 00_make, 12a, 02, 06, 22, 11 plus its helper, 24 and the Longcycler rebuild. Parsnp/Panaroo and script 14 are mentioned on their dedicated slides because they clarify population context or statistical limits. The comprehensive register remains an audit reference, not presentation content.",
  14: "Be explicit about assumptions instead of hiding them: the clinical rule is valid; one selected assembly represents the sampled episode; missing WGS and assembler fallback may be selective; adjacent samples are not continuous observation; and <=25 is only an operational support rule. Recover the culture-to-colony process, sequencing chemistry, basecalling/read QC/coverage, assembler and polishing commands, database releases, software versions, code revision and R environment. Keep antibiotic, demographic, host, wgMLST-reproduction and causal claims outside this evidence boundary.",
};

async function writeBlob(filePath, blob) {
  await fs.writeFile(filePath, new Uint8Array(await blob.arrayBuffer()));
}

function deleteShapeOrders(slide, orders) {
  const ids = orders
    .map((order) => slide.shapes.items[order - 1]?.id)
    .filter((id) => id !== undefined && id !== null);
  for (const id of ids.reverse()) slide.shapes.deleteById(id);
}

function addTextBox(slide, textValue, position, options = {}) {
  const shape = slide.shapes.add({ geometry: "textbox", position });
  shape.text.set(textValue);
  shape.text.style = {
    fontSize: options.fontSize ?? 14,
    color: options.color ?? "#111827",
    bold: options.bold ?? false,
  };
  shape.text.alignment = options.align ?? "left";
  shape.text.verticalAlignment = options.verticalAlign ?? "middle";
  return shape;
}

function addCard(slide, position, options = {}) {
  return slide.shapes.add({
    geometry: "roundRect",
    position,
    fill: options.fill ?? "#EEF3F8",
    line: {
      style: options.lineStyle ?? "solid",
      fill: options.stroke ?? "#CBD5E1",
      width: options.lineWidth ?? 1,
    },
  });
}

function connectBehind(slide, source, target, options = {}) {
  const connector = slide.shapes.connect(source, target, {
    kind: options.kind ?? "elbow",
    fromSide: options.fromSide ?? "right",
    toSide: options.toSide ?? "left",
    line: {
      style: options.dashed ? "dashed" : "solid",
      fill: options.color ?? "#94A3B8",
      width: options.width ?? 1.7,
    },
    tail: { type: "triangle", width: "sm", length: "sm" },
  });
  slide.elements.sendToBack(connector);
  return connector;
}

function addTopFlowNode(slide, x, count, label, caption, palette = {}) {
  const background = addCard(slide, { left: x, top: 175, width: 210, height: 92 }, {
    fill: palette.fill ?? "#E8F1FA",
    stroke: palette.stroke ?? "#AFC9E6",
  });
  addTextBox(slide, count, { left: x + 12, top: 184, width: 186, height: 28 }, {
    fontSize: 25, bold: true, color: palette.accent ?? "#2B6CB0", align: "center",
  });
  addTextBox(slide, label, { left: x + 12, top: 214, width: 186, height: 22 }, {
    fontSize: 13.5, bold: true, align: "center",
  });
  addTextBox(slide, caption, { left: x + 10, top: 238, width: 190, height: 20 }, {
    fontSize: 9.5, color: "#5F718D", align: "center",
  });
  return background;
}

function rebuildMethodsFlowchart(slide) {
  // Reuse the nine inherited flow cards in their original template positions.
  // Only the connector logic changes: the clinical track joins the assembly
  // track at the selected episode-level dataset, rather than appearing as one
  // sequential attrition funnel.
  const groupOrders = [
    [10, 11, 12, 13, 14, 15],
    [17, 18, 19, 20, 21, 22],
    [24, 25, 26, 27, 28, 29],
    [30, 31, 32, 33, 34, 35],
    [37, 38, 39, 40, 41, 42],
    [44, 45, 46, 47, 48, 49],
    [50, 51, 52, 53, 54, 55],
    [57, 58, 59, 60, 61, 62],
    [64, 65, 66, 67, 68, 69],
  ];
  const groups = groupOrders.map((orders) => orders.map((order) => slide.shapes.items[order - 1]));
  if (groups.flat().some((shape) => !shape)) throw new Error("Slide 4 inherited flow-card structure is incomplete");

  const inheritedArrows = [16, 23, 36, 43, 56, 63]
    .map((order) => slide.shapes.items[order - 1])
    .filter(Boolean);
  for (const arrow of inheritedArrows.reverse()) slide.shapes.deleteById(arrow.id);

  function rewriteCard(parts, copy) {
    const [background, title, detail, badge, badgeText, stage] = parts;
    title.text.set(copy.title);
    detail.text.set(copy.detail);
    badgeText.text.set(copy.badge);
    stage.text.set(copy.stage);
    return background;
  }

  addTextBox(slide, "ASSEMBLY / QC TRACK", { left: 70, top: 146, width: 260, height: 18 }, {
    fontSize: 11.5, bold: true, color: "#2B6CB0",
  });

  const copies = [
    { title: "1,303 FASTAs", detail: "Candidate files discovered on disk", badge: "1", stage: "discovery manifest" },
    { title: "1,299 metadata rows", detail: "1,295 matched + 4 audit-only; 8 discovered FASTAs absent", badge: "2", stage: "00_make reconciliation" },
    { title: "1,291 candidates", detail: "Eight rows outside primary/genomics curation; 579 episode keys", badge: "3", stage: "Longcycler/Flye alternatives" },
    { title: "1,211 QC pass", detail: "Eighty fail implemented QC; 556 episode keys remain", badge: "4", stage: "12a_wgs_qc.R" },
    { title: "556 selected genomes", detail: "One QC-pass genome per key; 655 passing alternatives not selected", badge: "5", stage: "canonical selection" },
    { title: "Clinically linked set", detail: "556 rows · 162 participants · 17 UTI · 539 Not_UTI; 532 Longcycler + 24 Flye", badge: "6", stage: "22 joins participant + timepoint" },
    { title: "Independent clinical track", detail: "583 episodes · 166 participants · 18 UTI · 565 Not_UTI; 27 lack selected WGS (4 + 23)", badge: "7", stage: "00a + 00b · UTI_Status" },
    { title: "Executed mixed route", detail: "Keep all 556 selected genomes; build 394 adjacent transitions", badge: "8", stage: "main analysis" },
    { title: "Longcycler sensitivity", detail: "Retain 532; exclude 24 Flye; rebuild 371 transitions with 9 new adjacencies", badge: "9", stage: "separate sensitivity route" },
  ];
  const cards = copies.map((copy, index) => rewriteCard(groups[index], copy));

  connectBehind(slide, cards[0], cards[1], { kind: "straight", width: 1.5 });
  connectBehind(slide, cards[1], cards[2], { kind: "straight", width: 1.5 });
  connectBehind(slide, cards[2], cards[3], { fromSide: "bottom", toSide: "top", kind: "straight", width: 1.5 });
  connectBehind(slide, cards[3], cards[4], { kind: "straight", width: 1.5 });
  connectBehind(slide, cards[4], cards[5], { kind: "straight", width: 1.5 });

  // The orange clinical connector enters from a separate track. Its shallow
  // diagonal occupies the gap between the middle and bottom card rows.
  connectBehind(slide, cards[6], cards[5], {
    fromSide: "top", toSide: "bottom", kind: "straight", color: "#D97706", width: 1.7,
  });
  connectBehind(slide, cards[5], cards[7], {
    fromSide: "bottom", toSide: "top", color: "#2B6CB0", width: 1.5,
  });
  connectBehind(slide, cards[5], cards[8], {
    fromSide: "bottom", toSide: "top", dashed: true, color: "#D97706", width: 1.5,
  });

}

function addParallelBranchSource(slide) {
  deleteShapeOrders(slide, [19, 20]);
  const source = addCard(slide, { left: 520, top: 148, width: 240, height: 28 }, {
    fill: "#EEF3F8", stroke: "#D7E0EA",
  });
  source.text.set("556 selected canonical genomes");
  source.text.style = { fontSize: 12.5, color: "#334155", bold: true };
  source.text.alignment = "center";
  source.text.verticalAlignment = "middle";
  const cards = [10, 13, 16].map((order) => slide.shapes.items[order - 1]).filter(Boolean);
  for (const card of cards) {
    connectBehind(slide, source, card, { fromSide: "bottom", toSide: "top", width: 1.3 });
  }
}

function tidyDistanceMethodSlide(slide) {
  const chipPairs = [
    [82, 83, 92, 112], [84, 85, 220, 112],
    [86, 87, 514, 112], [88, 89, 642, 112],
    [90, 91, 936, 112], [92, 93, 1060, 140],
  ];
  const chipObjects = chipPairs.map(([backgroundOrder, textOrder, left, width]) => ({
    background: slide.shapes.items[backgroundOrder - 1],
    textShape: slide.shapes.items[textOrder - 1],
    left,
    width,
  }));
  deleteShapeOrders(slide, [19, 20]);
  for (const { background, textShape, left, width } of chipObjects) {
    if (!background || !textShape) continue;
    background.position = { left, top: 336, width, height: 27 };
    textShape.position = { left: left + 5, top: 342, width: width - 10, height: 15 };
    if (width > 120) {
      textShape.text.style = { fontSize: 9, color: "#5F718D" };
      textShape.text.alignment = "center";
    }
  }
}

function tidyThresholdSlide(slide) {
  const title = slide.shapes.items[48];
  const body = slide.shapes.items[49];
  if (title) title.position = { left: 211.2, top: 596, width: 858.24, height: 22 };
  if (body) body.position = { left: 211.2, top: 618, width: 858.24, height: 20 };
}

await fs.mkdir(previewDir, { recursive: true });
await fs.mkdir(finalLayoutDir, { recursive: true });
await fs.mkdir(path.dirname(finalPptx), { recursive: true });

const presentation = await PresentationFile.importPptx(await FileBlob.load(starter));

for (let i = 1; i <= presentation.slides.items.length; i += 1) {
  const slide = presentation.slides.items[i - 1];
  const layoutPath = path.join(layoutDir, `starter-slide-${String(i).padStart(2, "0")}.layout.json`);
  const layout = JSON.parse(await fs.readFile(layoutPath, "utf8"));
  const textElements = layout.elements.filter((element) => typeof element.text === "string");
  const replacements = content[i];
  if (!replacements || replacements.length !== textElements.length) {
    throw new Error(`Slide ${i}: expected ${textElements.length} replacement strings, received ${replacements?.length ?? 0}`);
  }
  for (let j = 0; j < textElements.length; j += 1) {
    const element = textElements[j];
    const target = slide.shapes.items[element.order - 1];
    if (!target?.text) throw new Error(`Slide ${i}: cannot resolve source shape order ${element.order}`);
    const currentText = target.text.toString();
    if (currentText !== element.text) {
      throw new Error(`Slide ${i}: shape order ${element.order} contains '${currentText.slice(0, 60)}', expected '${element.text.slice(0, 60)}'`);
    }
    target.text.set(replacements[j]);
    const cardTitleOrders = new Set([11, 18, 25, 31, 38, 45, 51, 58, 65]);
    const shift = ((i === 4 || i === 13) && cardTitleOrders.has(element.order))
      ? 38
      : (i > 1 && element.order === 8)
        ? 22
        : 0;
    if (shift > 0) {
      const pos = target.position;
      target.position = {
        left: pos.left + shift,
        top: pos.top,
        width: Math.max(20, pos.width - shift),
        height: pos.height,
      };
    }
  }
  if (i === 1) deleteShapeOrders(slide, Array.from({ length: 12 }, (_, index) => index + 18));
  if (i === 4) rebuildMethodsFlowchart(slide);
  if (i === 6) addParallelBranchSource(slide);
  if (i === 10) tidyThresholdSlide(slide);
  if (i === 11) tidyDistanceMethodSlide(slide);
  slide.speakerNotes.textFrame.setText(notes[i]);
  slide.speakerNotes.setVisible(true);

  const preview = await presentation.export({ slide, format: "png", scale: 1.5 });
  await writeBlob(path.join(previewDir, `slide-${String(i).padStart(2, "0")}.png`), preview);
  const finalLayout = await slide.export({ format: "layout" });
  await fs.writeFile(path.join(finalLayoutDir, `slide-${String(i).padStart(2, "0")}.layout.json`), await finalLayout.text());
}

const montage = await presentation.export({ format: "webp", montage: true, scale: 1 });
await writeBlob(path.join(tmp, "final-montage.webp"), montage);
const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(finalPptx);

const inspected = await presentation.inspect({
  kind: "slide,textbox,shape,notes",
  maxChars: 200000,
});
await fs.writeFile(path.join(tmp, "final-inspect.ndjson"), inspected.ndjson);
console.log(finalPptx);

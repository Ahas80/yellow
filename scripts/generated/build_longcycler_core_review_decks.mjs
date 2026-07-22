#!/usr/bin/env node

import crypto from "node:crypto";
import { execFile } from "node:child_process";
import fs from "node:fs/promises";
import path from "node:path";
import { promisify } from "node:util";

import {
  assertNoEmptyPlaceholders,
  assertPptxHasNoForbiddenContent,
  bootstrapToScratch,
  loadRegistry,
  parseArgs,
  releaseCounts,
  releaseDate,
  writeCsvAtomic,
  writeJsonAtomic,
} from "./longcycler_release_presentation_common.mjs";

const execFileAsync = promisify(execFile);
const SKILL_DIR = "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.709.11516/skills/presentations";
const PYTHON = "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";
const DEFAULT_PROJECT_ROOT = "/Users/Aamir/Desktop/rUTIs";
const RETIRED_INPUT_TOKEN = ["fl", "ye"].join("");
const TRANSPARENT_PNG = Buffer.from(
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M/wHwAF/gL+Wf6ZVQAAAABJRU5ErkJggg==",
  "base64",
);

const DECKS = [
  {
    key: "onboarding",
    label: "clinical-genomic onboarding",
    source: "outputs/manual-20260526-ruti-onboarding/presentations/ruti-clinical-genomic-onboarding/output/ruti-clinical-genomic-onboarding-review.pptx",
    output: "outputs/manual-20260526-ruti-onboarding/presentations/ruti-clinical-genomic-onboarding/output/ruti-clinical-genomic-onboarding-review.pptx",
    slides: 15,
    imagePlan: {
      6: "plots/vf/uti_not_uti_denominator_waterfall.png",
      8: "plots/mechanism/not_uti_to_uti_case_matrix.png",
      9: "plots/mechanism/host_context_transition_heatmap.png",
      10: "plots/final/Fig05_vf_association_evidence.png",
      13: "plots/final/supplementary/FigS07_plasmid_amr_context.png",
      14: "plots/final/Fig08_reference_aware_variant_map.png",
      15: "plots/final/Fig06_longitudinal_trajectories.png",
    },
  },
  {
    key: "scientific",
    label: "scientific review",
    source: "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review/output/rUTI_Current_Results_Scientific_Review_2026-05-27.pptx",
    output: "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review/output/rUTI_Current_Results_Scientific_Review_2026-05-27.pptx",
    slides: 22,
    nativeSlides: new Set([6, 14]),
    imagePlan: {
      4: "plots/vf/uti_not_uti_denominator_waterfall.png",
      5: "plots/clinical/uti_not_uti_clinical_rule_flow.png",
      7: "plots/mechanism/not_uti_to_uti_case_matrix.png",
      8: "plots/mechanism/not_uti_to_uti_case_matrix.png",
      9: "plots/mechanism/host_context_transition_heatmap.png",
      10: "plots/final/Fig05_vf_association_evidence.png",
      15: "plots/clinical/uti_not_uti_clinical_rule_flow.png",
      16: "plots/final/supplementary/FigS06_transition_mechanisms.png",
      17: "plots/final/supplementary/FigS07_plasmid_amr_context.png",
      18: "plots/final/supplementary/FigS05_near_miss_leave_one_uti.png",
      19: "plots/vf/uti_not_uti_leave_one_uti_out_stability.png",
      20: "plots/final/Fig08_reference_aware_variant_map.png",
      21: "plots/final/Fig06_longitudinal_trajectories.png",
    },
  },
  {
    key: "vf_focused",
    label: "VF-focused review",
    source: "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx",
    output: "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx",
    slides: 22,
    imagePlan: {
      4: "plots/vf/uti_not_uti_denominator_waterfall.png",
      7: "plots/vf/vf_burden_boxplot.png",
      8: "plots/vf/vf_gene_screening_vs_model_evidence.png",
      10: "plots/mechanism/not_uti_to_uti_case_matrix.png",
      11: "plots/mechanism/host_context_transition_heatmap.png",
      12: "plots/final/Fig05_vf_association_evidence.png",
      15: "plots/clinical/uti_not_uti_clinical_rule_flow.png",
      16: "plots/vf/module_prevalence_by_status.png",
      17: "plots/final/supplementary/FigS06_transition_mechanisms.png",
      18: "plots/final/supplementary/FigS07_plasmid_amr_context.png",
      19: "plots/final/supplementary/FigS05_near_miss_leave_one_uti.png",
      20: "plots/vf/uti_not_uti_leave_one_uti_out_stability.png",
      21: "plots/final/Fig08_reference_aware_variant_map.png",
    },
  },
];

async function runTool(command, args, options = {}) {
  const { stdout, stderr } = await execFileAsync(command, args, {
    cwd: options.cwd,
    env: { ...process.env, ...(options.env || {}) },
    encoding: "utf8",
    maxBuffer: 128 * 1024 * 1024,
  });
  if (options.echo && stdout.trim()) console.log(stdout.trim());
  if (options.echo && stderr.trim()) console.error(stderr.trim());
  return { stdout, stderr };
}

function parseOverflowSlides(output) {
  const match = String(output || "").match(/Slides with content overflowing original canvas \(1-based indexing\):\s*([^\n]+)/i);
  if (!match) return [];
  return match[1]
    .split(",")
    .map((value) => Number.parseInt(value.trim(), 10))
    .filter(Number.isInteger);
}

async function auditOverflowAgainstStarter({ deck, workspace, starterPptx, finalPptx }) {
  const testScript = path.join(SKILL_DIR, "container_tools", "slides_test.py");
  const starterResult = await runTool(PYTHON, [testScript, starterPptx]);
  const finalResult = await runTool(PYTHON, [testScript, finalPptx]);
  const starterSlides = parseOverflowSlides(`${starterResult.stdout}\n${starterResult.stderr}`);
  const finalSlides = parseOverflowSlides(`${finalResult.stdout}\n${finalResult.stderr}`);
  const starterSet = new Set(starterSlides);
  const finalSet = new Set(finalSlides);
  const introducedSlides = finalSlides.filter((slide) => !starterSet.has(slide));
  const resolvedSlides = starterSlides.filter((slide) => !finalSet.has(slide));
  const audit = {
    deck: deck.key,
    status: introducedSlides.length ? "FAIL" : "PASS",
    starter_overflow_slides: starterSlides,
    final_overflow_slides: finalSlides,
    introduced_overflow_slides: introducedSlides,
    resolved_overflow_slides: resolvedSlides,
    interpretation: starterSlides.length
      ? "Inherited template exceptions are accepted only when the final deck introduces no additional overflow slides."
      : "The starter and final deck must both pass with no overflow slides.",
  };
  await writeJsonAtomic(path.join(workspace, "overflow-audit.json"), audit);
  if (introducedSlides.length) {
    throw new Error(`${deck.key}: final deck introduced overflow on slides ${introducedSlides.join(", ")}`);
  }
  if (!starterSlides.length && finalSlides.length) {
    throw new Error(`${deck.key}: final deck failed overflow audit on slides ${finalSlides.join(", ")}`);
  }
  console.log(
    `${deck.key}: overflow audit PASS; starter=${starterSlides.length}, final=${finalSlides.length}, introduced=0`,
  );
}

function slidesFromPresentation(presentation) {
  if (Array.isArray(presentation.slides?.items)) return presentation.slides.items;
  if (Number.isInteger(presentation.slides?.count) && typeof presentation.slides.getItem === "function") {
    return Array.from({ length: presentation.slides.count }, (_, index) => presentation.slides.getItem(index));
  }
  throw new Error("Could not enumerate imported slides.");
}

async function writeBlob(filePath, blob) {
  await fs.mkdir(path.dirname(filePath), { recursive: true });
  await fs.writeFile(filePath, new Uint8Array(await blob.arrayBuffer()));
}

function fileBytes(buffer) {
  return buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength);
}

function hashBuffer(buffer) {
  return crypto.createHash("sha256").update(buffer).digest("hex");
}

function withReleaseDate(text, date) {
  return String(text)
    .replaceAll("27 May 2026", date)
    .replaceAll("26 May 2026", date)
    .replaceAll("current UTI_Status analysis", "operational UTI phenotype | Longcycler-only analysis");
}

function lookup(map, slide, source, date, aid, occurrenceIndex = 0) {
  const replacement =
    map[slide]?.__byAid?.[aid] ??
    map[slide]?.__byOccurrence?.[source]?.[occurrenceIndex] ??
    map[slide]?.[source];
  return withReleaseDate(replacement === undefined ? source : replacement, date);
}

function setInheritedTextFrame(target, bbox, { width = bbox?.[2], height = bbox?.[3] } = {}) {
  if (!bbox || bbox.length !== 4) return;
  target.position = {
    left: bbox[0],
    top: bbox[1],
    width,
    height,
  };
}

function applyInheritedTextVisualRepair({ deckKey, slideNumber, target, element, next }) {
  if (deckKey !== "onboarding" || !target?.text) return next;

  let repairedText = next;
  if (slideNumber === 11 && element.bbox?.[1] >= 530 && element.bbox?.[1] < 536 && element.bbox?.[2] < 30) {
    repairedText = ">";
    if (repairedText !== next) target.text.set(repairedText);
  } else if (slideNumber === 11 && element.bbox?.[1] >= 536 && element.bbox?.[1] <= 545 && element.bbox?.[2] >= 100) {
    if (element.bbox[0] >= 700 && element.bbox[0] < 900) {
      repairedText = "scripts/research_questions/run_all.R";
    } else if (element.bbox[0] >= 900) {
      repairedText = "RUN_COMPLETE_ANALYSIS.sh";
    }
    if (repairedText !== next) target.text.set(repairedText);
  }

  if (slideNumber === 4 && next === "Operational UTI requires culture and compatible symptoms.") {
    target.text.fontSize = 32;
  }
  if (slideNumber === 5 && next === "not_uti_to_uti_casebook.csv") {
    target.text.fontSize = 8;
    setInheritedTextFrame(target, element.bbox, { width: 190 });
  }
  if (slideNumber === 5 && next === "RUN_COMPLETE_ANALYSIS.sh") {
    target.text.fontSize = 8;
    setInheritedTextFrame(target, element.bbox, { width: 220 });
  }
  if (slideNumber === 6 && next === "16") {
    target.text.fontSize = 43;
    setInheritedTextFrame(target, element.bbox, { width: 100, height: 58 });
  }
  if (slideNumber === 6 && next === "516") {
    setInheritedTextFrame(target, element.bbox, { width: 100 });
  }
  if (slideNumber === 7 && next === "0 endpoints missing") {
    setInheritedTextFrame(target, element.bbox, { width: 160 });
  }
  if (slideNumber === 10 && next === "16") {
    target.text.fontSize = 43;
    setInheritedTextFrame(target, element.bbox, { width: 100, height: 58 });
  }
  if (slideNumber === 10 && next === "532") {
    setInheritedTextFrame(target, element.bbox, { width: 100 });
  }
  return repairedText;
}

function onboardingTextMap(c) {
  return {
    1: {
      "SCIENTIFIC ONBOARDING REVIEW": "LONGCYCLER-ONLY SCIENTIFIC ONBOARDING",
      "Recurrent UTI clinical-to-genomic analysis in E. coli": "Longcycler-only clinical-to-genomic analysis of the operational UTI phenotype",
      "CURRENT EVIDENCE": "CURRENT RELEASE",
      "583": String(c.episodes),
      "Clinical episodes": "Selected episodes",
      "Current primary denominator": `${c.residents} residents`,
      "18": String(c.uti),
      "Primary UTI": "Operational UTI",
      "11": String(c.toUti),
      "Current framing: UTI versus heterogeneous Not_UTI.\nLegacy ASB comparisons are not primary conclusions.": "Current framing: operational UTI versus heterogeneous Not_UTI.\nAll analytical inputs are selected QC-passing Longcycler genomes.",
      "Current evidence pack: results/summary/final_key_results_summary.md": "Current evidence pack: results/pipeline/longcycler_release_claim_registry.json",
    },
    2: {
      "The primary contrast asks which episodes fulfil a clinical UTI rule, not whether E. coli is present.": "The primary contrast uses a versioned operational UTI rule; it does not equate bacterial presence with infection.",
      "18 current primary clinical episodes": `${c.uti} operational UTI episodes`,
      "565 episodes; clinically heterogeneous": `${c.notUti} selected episodes; clinically heterogeneous`,
      "Source: current primary UTI_Status definition; validation checks (18 UTI / 565 Not_UTI)": `Source: operational UTI phenotype; selected cohort (${c.uti} UTI / ${c.notUti} Not_UTI)`,
    },
    3: {
      "583": String(c.episodes),
      "Clinical episodes": "Selected episodes",
      "167 participants": `${c.residents} residents`,
      "556": String(c.episodes),
      "VF/WGS-ready": "QC-passing genomes",
      "Linked genomic episodes": "one selected per episode",
      "227": String(c.vfFeatures),
      "VF gene columns": "VFDB features",
      "11": String(c.transitions),
      "Transitions": "Adjacent pairs",
      "Not_UTI to UTI": `${c.transitionResidents} residents; ${c.toUti} to UTI`,
      "Sources: final_key_results_summary.md; final_figure_validation_checks.csv": "Source: results/pipeline/longcycler_release_claim_registry.json",
    },
    4: {
      "A UTI episode must satisfy two clinical conditions.": "Operational UTI requires culture and compatible symptoms.",
      "The operational UTI phenotype requires two clinical conditions.": "Operational UTI requires culture and compatible symptoms.",
      "The operational UTI phenotype requires culture support and compatible symptoms.": "Operational UTI requires culture and compatible symptoms.",
      "Operational UTI requires culture support and compatible symptoms.": "Operational UTI requires culture and compatible symptoms.",
      "v": "↓",
      "583": String(c.episodes),
      "primary clinical episodes assessed": "selected Longcycler-linked episodes assessed",
      "UTI  |  18 episodes": `UTI  |  ${c.uti} episodes`,
      "Not_UTI  |  565 episodes": `Not_UTI  |  ${c.notUti} episodes`,
      "This definition sets the evidence boundary for every genomic result that follows.": "This versioned operational definition sets the boundary for every genomic result that follows.",
      "Source: current primary rule; final_figure_validation_checks.csv (n=583, UTI=18, Not_UTI=565)": `Source: selected cohort (n=${c.episodes}, UTI=${c.uti}, Not_UTI=${c.notUti})`,
    },
    5: {
      "Input batches": "Operational phenotype",
      "00a_load_clean_clinical.R": "00a / 00b clinical classification",
      "Primary status": "Selected cohort",
      "00b_classify_episodes.R": "analysis_cohort_longcycler.csv",
      "status_map.csv": "Participant_id + tp_lab",
      "Assemblies": "Selected assemblies",
      "FASTA / GFF inputs": "QC-passing Longcycler only",
      "Genomic profile": "Genome evidence",
      "VF + ST + SNPs": "VF + ST + direct SNP pairs",
      "strain_compare": "pairwise_metrics.csv",
      "not_uti_to_uti_casebook.csv": "not_uti_to_uti_casebook.csv",
      "results/mechanism/not_uti_to_uti_casebook.csv": "not_uti_to_uti_casebook.csv",
      "Final figure pack: denominator -> transitions -> host context -> robustness": "Release layer: denominator -> direct pairs -> transitions -> RQ01-RQ10",
      "35_final_figure_pack.R": "RUN_COMPLETE_ANALYSIS.sh",
      "Source: docs/pipeline_architecture.md; canonical outputs and scripts named on slide": "Source: central selected cohort, claim registry and final-run evidence products",
    },
    6: {
      "The dataset is rich; the UTI denominator is small.": "Rich dataset; sparse operational UTI outcome.",
      "17": String(c.uti),
      "UTI episodes in the\nVF/model-ready cohort": "operational UTI episodes\nin the selected cohort",
      "539": String(c.notUti),
      "Source: Main Figure 1; final_figure_validation_checks.csv (validated denominators)": "Source: current final-run denominator figure; selected cohort only",
    },
    7: {
      __byOccurrence: {
        "1.00": [String(c.directPairs), String(c.transitions)],
        "0 / 0": [String(c.transitionsAtThreshold), String(c.casebookLinked)],
      },
      "Participant 20026 changes clinical state without measured profile acquisition.": `${c.toUti} operational Not_UTI-to-UTI transitions are all genomically linked.`,
      "One concrete Not_UTI to UTI transition anchors the aggregate interpretation.": `Direct pair evidence is primary: ${c.toUtiAtThreshold} of ${c.toUti} are at or below ${c.snpThreshold} SNPs.`,
      "DESCRIPTIVE WORKED EXAMPLE": "DESCRIPTIVE AGGREGATE",
      "T3": `${c.toUti} cases`,
      "02 Dec 2024": `${c.transitions} adjacent pairs`,
      "10,000-100,000 CFU/mL": `${c.transitionResidents} residents`,
      "No defining symptoms recorded": "Operational Not_UTI endpoint",
      "42 DAYS": `${c.snpThreshold} SNP CUTOFF`,
      "5 SNPS | STRONG LINK": `${c.toUtiAtThreshold} OF ${c.toUti} AT OR BELOW`,
      "UTI-1": `${c.casebookLinked} linked`,
      "13 Jan 2025": `${c.casebookCases} casebook rows`,
      ">100,000 CFU/mL": `${c.casebookMissing} missing endpoints`,
      "0 missing endpoints": "0 endpoints missing",
      "Symptoms now satisfy UTI definition": "Operational UTI endpoint",
      "VF JACCARD": "ALL DIRECT PAIRS",
      "MODULE JACCARD": "ADJACENT PAIRS",
      "PLASMID GAIN / LOSS": `PAIRS <= ${c.snpThreshold} SNPs`,
      "AMR GAIN / LOSS": "CASEBOOK LINKED",
      "Symptoms": "Operational",
      "CLINICAL CHANGE": "PHENOTYPE",
      "Stable measured bacterial profile + changing clinical state: host context or unmeasured expression/regulation remains plausible.": "Close direct pairs motivate host-state or regulatory hypotheses; distant pairs remain replacement-consistent. Neither pattern proves mechanism.",
      "Source: not_uti_to_uti_casebook.csv; participant 20026 / case_191": "Source: direct pair table, canonical transitions and 9-row linked casebook",
    },
    8: {
      "Eleven clinical transitions are tracked; ten have linked WGS/VF evidence.": `${c.toUti} operational transitions are tracked; all ${c.casebookLinked} have linked genomic evidence.`,
      "10 WGS/VF LINKED": `${c.casebookLinked}/${c.casebookCases} GENOMICALLY LINKED`,
      "4  same-strain,\nstable-profile": `${c.toUtiAtThreshold}  at or below\n${c.snpThreshold} SNPs`,
      "3  consistent with\nstrain replacement": `${c.toUti - c.toUtiAtThreshold}  above the\noperational threshold`,
      "2  uncertain\nmechanisms": `${c.directPairs} direct pairs\nin the selected cohort`,
      "1 profile-change case\n1 missing genomic endpoint": `${c.casebookMissing} missing genomic endpoints\nmechanism labels are descriptive`,
      "Source: Main Figure 2; final_figure_validation_checks.csv and captions": "Source: current final-run mechanism casebook and direct pair evidence",
    },
    9: {
      "Stable bacterial profiles can coexist with changing clinical state.": "Direct paired evidence separates close and distant clinical transitions.",
      "Low SNP distance and stable VF content do not guarantee the same symptom state.": `${c.toUtiAtThreshold} of ${c.toUti} focused transitions are at or below ${c.snpThreshold} SNPs; ${c.toUti - c.toUtiAtThreshold} are above.`,
      "4": String(c.toUtiAtThreshold),
      "same-strain stable-profile\ntransitions": `focused transitions <=${c.snpThreshold} SNPs`,
      "Persistent strain evidence moves the next hypothesis toward host state or unmeasured regulation.": "Close direct pairs move the next hypothesis toward host state or unmeasured regulation.",
      "Stable profile does not prove causality, persistence, or a non-genomic trigger.": "A distance threshold is operational evidence, not mechanistic proof.",
      "Source: Main Figure 3; final_figure_captions.md (hypothesis-generating)": "Source: current final-run direct pair and host-context figure",
    },
    10: {
      "Unadjusted patterns are contextual only: no association remains significant after FDR correction.": `The ${c.uti}-episode operational UTI group supports exploratory modelling only; no causal claim is made.`,
      "0": String(c.uti),
      "FDR-significant global\nVF associations": "operational UTI episodes\nin the selected cohort",
      "17": String(c.episodes),
      "VF-ready UTI episodes": "selected Longcycler-linked episodes",
      "Do not revive older ASB/UTI or significant lpf narratives as current conclusions.": "Do not promote exploratory VF screens, named mutations or AMR context to causal conclusions.",
      "Source: Main Figure 4; final_figure_captions.md (exploratory; sparse UTI denominator)": "Source: current final-run VF robustness figure; operational phenotype",
    },
    11: {
      __byOccurrence: {
        "35_final_figure_pack.R": ["scripts/research_questions/run_all.R", "RUN_COMPLETE_ANALYSIS.sh"],
        "scripts/research_questions/run_all.R": [
          "scripts/research_questions/run_all.R",
          "scripts/research_questions/run_all.R",
          "RUN_COMPLETE_ANALYSIS.sh",
        ],
      },
      "Five output families locate the evidence; five scripts regenerate it.": "One claim registry locates the evidence; the complete runner regenerates RQ01-RQ10.",
      "A new colleague can start with the interpretation, then trace every figure to its data product.": `${c.episodes} episodes | ${c.directPairs} direct pairs | ${c.transitions} adjacent pairs from ${c.transitionResidents} residents | ${c.nearMiss} near-miss rows`,
      "results/summary/final_key_results_summary.md": "results/pipeline/longcycler_release_claim_registry.json",
      "results/final_figures/final_figure_captions.md": "results/research_questions/RUN_COMPLETE.txt",
      "Figure captions + limits": "RQ release marker",
      "results/final_figures/final_figure_manifest.csv": "results/strain_compare/pairwise_metrics.csv",
      "Validated figure inventory": "Direct-pair evidence",
      "results/mechanism/not_uti_to_uti_casebook.md": "results/mechanism/not_uti_to_uti_casebook.csv",
      "results/clinical/status_map.csv": "results/clinical/analysis_cohort_longcycler.csv",
      "results/vf/vf_analysis_ready.csv": "results/qc/analysis_assembly_manifest.csv",
      "VF analysis cohort": "Genome manifest",
      "Selected genome manifest": "Genome manifest",
      "docs/pipeline_architecture.md": "RUN_COMPLETE_ANALYSIS.sh",
      "Pipeline explanation": "Complete runner",
      "Complete pipeline runner": "Complete runner",
      "35_final_figure_pack.R": "scripts/research_questions/run_all.R",
      "Final figure script": "RQ01-RQ10 runner",
      "CLASSIFY": "SELECT COHORT",
      "00b_classify_episodes.R": "analysis_cohort_longcycler.csv",
      "BUILD VF COHORT": "DIRECT EVIDENCE",
      "22_vf_build_analysis_dataset.R": "pairwise_metrics.csv",
      "MECHANISMS": "TRANSITIONS",
      "33_mechanism_first_addon.R": "longcycler_transitions.csv",
      "ROBUSTNESS": "RQ01-RQ10",
      "34_robustness_first_addon.R": "scripts/research_questions/run_all.R",
      "FIGURES": "RELEASE",
      "Canonical paths verified in workspace on 26 May 2026": "Canonical release paths are verified before deck generation",
    },
    12: {
      "The next question is what triggers symptoms in a persisting strain.": "The next question is what distinguishes close direct pairs from clinical symptom emergence.",
      "4 stable-profile transitions": `${c.toUtiAtThreshold} focused pairs <=${c.snpThreshold} SNPs`,
      "Same-strain evidence makes symptom emergence without measured acquisition a real, testable pattern.": "Close direct pairs make host-state and regulation hypotheses testable without proving mechanism.",
      "17 UTI genomic episodes": `${c.uti} operational UTI episodes`,
      "Which transition deserves deeper follow-up?": `Which of the ${c.toUti} linked transitions deserves deeper follow-up?`,
      "Primary story: UTI versus Not_UTI | 22 min review + 8 min discussion": "Primary story: operational UTI versus Not_UTI | Longcycler-only release",
    },
    13: {
      "Cached ResFinder screening supports the mechanism supplement at transition level.": "Current screening provides exploratory accessory, plasmid and AMR context for linked transitions.",
    },
    15: {
      "Source: results/longitudinal/swimmer_plot.png | contextual lookup only": "Source: plots/publication/Fig1_Swimmer_Plot.png | contextual lookup only",
    },
  };
}

function scientificTextMap(c) {
  return {
    1: {
      "rUTI Clinical-to-Genomic Analysis: Current UTI / Not_UTI Evidence": "rUTI Clinical-to-Genomic Analysis: Longcycler-Only Operational UTI Evidence",
      "Do strains persist? What changes when symptoms emerge? What can sparse UTI events support? | Scientific review and handover": `Selected cohort: ${c.episodes} episodes from ${c.residents} residents | RQ01-RQ10 scientific review and handover`,
    },
    2: {
      "Repeated bacteriuria does not automatically indicate symptomatic infection.\nUTI | Culture-supported episode with compatible, catheter-aware symptoms.\nNot_UTI | Comparator state; clinically heterogeneous by design.\nLongitudinal sampling asks what changes before symptoms emerge.\nPrimary analysis throughout: UTI versus Not_UTI.": "Repeated bacteriuria does not automatically indicate symptomatic infection.\nUTI | Versioned operational phenotype: culture support plus compatible, catheter-aware symptoms.\nNot_UTI | Comparator state; clinically heterogeneous by design.\nLongitudinal sampling asks what changes before symptoms emerge.\nPrimary analysis throughout: operational UTI versus Not_UTI in selected Longcycler-linked episodes.",
    },
    3: {
      "Repeated clinical episodes linked by participant and time.\nClinical lane: culture, symptoms and catheter context establish primary status.\nGenomic lane: isolate sequencing produces VF, SNP/ST and accessory context.\nIntegration asks whether clinical transitions reflect persistence, replacement or uncertainty.": `The selected cohort contains ${c.episodes} episodes from ${c.residents} residents.\nClinical lane: culture, symptoms and catheter context establish the operational phenotype.\nGenomic lane: one QC-passing Longcycler assembly per episode supplies VF, SNP/ST and accessory context.\nDirect pair evidence comprises ${c.directPairs} within-resident pairs and ${c.transitions} adjacent pairs.`,
    },
    4: {
      "DESCRIPTIVE | Clinical: 583 episodes (18 UTI; 565 Not_UTI). VF/model-ready: 556 (17 UTI; 539 Not_UTI).": `ANALYTICAL | Selected cohort: ${c.episodes} episodes (${c.uti} operational UTI; ${c.notUti} Not_UTI) from ${c.residents} residents.`,
    },
    5: {
      "Current UTI definition: culture support plus catheter-aware compatible symptoms": "Operational UTI phenotype: culture support plus catheter-aware compatible symptoms",
      "DESCRIPTIVE | Not_UTI is the comparator state, not a single biological phenotype.": `DESCRIPTIVE | Versioned operational rule applied to ${c.episodes} selected episodes; Not_UTI remains heterogeneous.`,
    },
    6: {
      "Integrated workflow: clinical status and genomic evidence meet at transition analysis": "Selected clinical status and direct genomic evidence meet at transition analysis",
      "WORKFLOW | Primary status, VF/WGS outputs and longitudinal casebook are integrated after episode linkage.": `WORKFLOW | ${c.episodes} selected episodes -> ${c.directPairs} direct pairs -> ${c.transitions} adjacent pairs -> RQ01-RQ10.`,
    },
    7: {
      "Worked transition: participant 20026 develops symptoms without measured profile change": `${c.toUti} operational Not_UTI-to-UTI transitions are all genomically linked`,
      "DESCRIPTIVE | Stable measured profile supports a host-state or regulatory hypothesis; it does not prove mechanism.": `DESCRIPTIVE | ${c.toUtiAtThreshold} of ${c.toUti} focused direct pairs are at or below ${c.snpThreshold} SNPs; no mechanism is proved.`,
    },
    8: {
      "DESCRIPTIVE | 11 clinical transitions; 10 WGS/VF-linked; 4 stable-profile and 3 replacement-consistent.": `DESCRIPTIVE | ${c.casebookCases} operational transitions; ${c.casebookLinked}/${c.casebookCases} genomically linked; ${c.casebookMissing} missing endpoint; ${c.toUtiAtThreshold} <=${c.snpThreshold} SNPs.`,
    },
    9: {
      "Central finding: the same measured strain profile can accompany changing clinical state": "Direct paired evidence distinguishes close from distant clinical transitions",
      "DESCRIPTIVE | Low SNP distance plus stable VF profile can accompany symptom emergence.": `DESCRIPTIVE | ${c.transitionsAtThreshold}/${c.transitions} adjacent pairs are <=${c.snpThreshold} SNPs; focused transitions: ${c.toUtiAtThreshold}/${c.toUti}.`,
    },
    10: {
      "Population-level VF findings remain exploratory under sparse UTI counts": `Population-level VF analysis remains exploratory with ${c.uti} operational UTI episodes`,
      "EXPLORATORY | No global VF association remains significant after FDR correction; 17 VF-ready UTI episodes.": "EXPLORATORY | Participant-aware models and robustness panels are descriptive and do not support causal inference.",
    },
    11: {
      "DESCRIPTIVE TRANSITION EVIDENCE\n11 clinical Not_UTI-to-UTI transitions; 10 have WGS/VF endpoints.\nFour stable-profile cases and three replacement-consistent cases shape the biological story.\n\nEXPLORATORY POPULATION EVIDENCE\nNo FDR-significant global VF association; sparse UTI counts constrain precision.\nSensitivity panels diagnose fragility rather than confirm association.": `DESCRIPTIVE DIRECT EVIDENCE\n${c.directPairs} within-resident pairs; ${c.transitions} adjacent pairs from ${c.transitionResidents} residents.\n${c.transitionsAtThreshold}/${c.transitions} adjacent pairs and ${c.toUtiAtThreshold}/${c.toUti} focused pairs are <=${c.snpThreshold} SNPs.\nCasebook: ${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing} cases/linked/missing.\n\nEXPLORATORY POPULATION EVIDENCE\n${c.uti} operational UTI episodes constrain precision.\n${c.nearMiss} near-miss rows remain sensitivity-only; RQ01-RQ10 reports guardrails.`,
    },
    12: {
      "WHAT THE DATA SUPPORT\nStable measured strain profiles can accompany symptom emergence.\nWHAT THEY DO NOT PROVE\nHost state, bacterial regulation/expression and unmeasured changes remain alternatives.\nPRIMARY LIMITATIONS\n18 clinical UTI episodes; 17 VF-ready UTI episodes; one transition missing genomic endpoint.\nRepeated measures and missing endpoints restrict population-level inference.": `WHAT THE DATA SUPPORT\nDirect genomic distance can identify close and distant clinical transitions.\nWHAT THEY DO NOT PROVE\nHost state, bacterial regulation/expression and unmeasured changes remain alternatives.\nPRIMARY LIMITATIONS\n${c.uti} operational UTI episodes; observational design; sparse outcome.\nAll ${c.casebookCases} focused transitions are linked, but repeated measures restrict population-level inference.`,
    },
    13: {
      "EVIDENCE | Preserve the transition-focused biological narrative.\nNEXT ANALYSIS | Prioritise within-strain follow-up and clinically grounded host/context variables.\nDESIGN | Analyse replacement separately from possible symptom emergence within strain.\nREPORTING | Keep global VF analysis exploratory; keep AMR as transition context only.\n\nDISCUSSION | What should become the next confirmatory question?": "EVIDENCE | Preserve the direct-pair and transition-focused narrative.\nNEXT ANALYSIS | Prioritise close within-strain pairs and clinically grounded host/context variables.\nDESIGN | Analyse distant replacement-consistent pairs separately.\nREPORTING | Keep VF and AMR results exploratory; use RQ01-RQ10 as the release index.\n\nDISCUSSION | What should become the next confirmatory question?",
    },
    22: {
      "CURRENT SOURCES OF TRUTH\nfinal_key_results_summary.md | headline results\nfinal_figure_captions.md + final_figure_manifest.csv | figures and limitations\nnot_uti_to_uti_casebook.md | transition case details\npipeline_architecture.md + 35_final_figure_pack.R | orientation and figure generation\nDECK UPDATE | Superseded result frames replaced with current UTI / Not_UTI evidence in an intentional academic-report restyle.\nGUARDRAILS | UTI versus Not_UTI only; no FDR-significant global VF claim; AMR is context only.\nDenominators displayed in this deck are validated episode-level counts.": `CURRENT SOURCES OF TRUTH\nlongcycler_release_claim_registry.json | validated claims and plot inventory\nanalysis_cohort_longcycler.csv | ${c.episodes} selected episodes\npairwise_metrics.csv | ${c.directPairs} direct pairs\nlongcycler_transitions.csv | ${c.transitions} adjacent pairs from ${c.transitionResidents} residents\nnot_uti_to_uti_casebook.csv | ${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing} cases/linked/missing\nRUN_COMPLETE.txt | RQ01-RQ10 complete\nGUARDRAILS | Operational UTI phenotype; Longcycler-only analytics; VF and AMR are exploratory.`,
    },
  };
}

function vfTextMap(c) {
  return {
    1: {
      "CURRENT RESULTS REVIEW  /  ONBOARDING": "LONGCYCLER-ONLY RESULTS REVIEW / ONBOARDING",
      "A current UTI versus Not_UTI review for a new project colleague": "Operational UTI versus Not_UTI evidence for a new project colleague",
    },
    2: {
      "The primary comparison separates culture-supported symptomatic episodes from a heterogeneous remainder.": "The versioned operational phenotype separates culture-supported compatible episodes from a heterogeneous remainder.",
    },
    3: {
      "The same sampling structure supports both clinical classification and VF/WGS comparisons.": `One selected QC-passing Longcycler genome links to each of ${c.episodes} analytical episodes.`,
      "Link episode status to WGS and virulence-factor profiles": "Link operational episode status to selected genomes and VF profiles",
    },
    4: {
      "The primary UTI denominator is deliberately narrow": "The selected operational UTI denominator is deliberately narrow",
    },
    5: {
      "Counts are shown at their analytical level: assembly-level QC is not an episode denominator.": "Every analytical count refers to the selected Longcycler-linked cohort; source attrition is excluded from results.",
      "Episodes classified": "Selected episodes",
      "585": String(c.episodes),
      "before primary manual exclusions": `${c.residents} residents`,
      "00a / 00b classification": "analysis_cohort_longcycler.csv",
      "Primary included": "Operational phenotype",
      "583": String(c.episodes),
      "18 UTI  |  565 Not_UTI  |  2 excluded": `${c.uti} UTI  |  ${c.notUti} Not_UTI`,
      "00b status_map_primary_included": "versioned operational rule",
      "Assembly QC records": "Selected assemblies",
      "1,291": String(c.episodes),
      "includes assembler alternatives": "QC-passing Longcycler only",
      "12a_wgs_qc.R": "analysis_assembly_manifest.csv",
      "Canonical episodes": "Selected genomes",
      "556": String(c.episodes),
      "one selected assembly per episode": "one genome per selected episode",
      "VF presence/absence": "VF feature rows",
      "episode-level VF rows": "selected episode-level rows",
      "Model-ready": "Direct-pair ready",
      "17 UTI  |  539 Not_UTI\n27 clinical episodes lack VF-ready evidence": `${c.directPairs} within-resident pairs\n${c.transitions} adjacent pairs`,
      "22_vf_build_analysis_dataset.R": "pairwise_metrics.csv",
      "VF feature space": "VFDB feature space",
      "227 -> 32": String(c.vfFeatures),
      "genes -> modules\n18 UPEC-candidate modules": "binary VF features\nmodule summaries exploratory",
      "26_vf_define_gene_modules.R": "current final-run VF outputs",
      "11 -> 10": `${c.toUti} -> ${c.casebookLinked}`,
      "clinical transitions -> WGS/VF-linked\n1 endpoint missing": `operational transitions -> genomically linked\n${c.casebookMissing} endpoints missing`,
    },
    6: {
      "556": String(c.episodes),
      "VF/WGS-linked\nepisodes": "selected Longcycler-linked\nepisodes",
      "227 VF genes": `${c.vfFeatures} VF features`,
      "32 modules": "VF modules",
      "18 UPEC-candidate systems": "exploratory summaries",
    },
    7: {
      "Observed VF burden overlaps strongly across status groups": "Observed VF burden is descriptive across operational status groups",
      "The UTI group is sparse (n=17), so this is descriptive distributional evidence rather than a causal test.": `The operational UTI group is sparse (n=${c.uti}), so this remains descriptive rather than causal.`,
    },
    8: {
      "Nominal VF screens do not become confirmed associations": "Participant-aware VF screens remain exploratory",
      "Participant-aware modelling and FDR correction define the evidential limit of the current dataset.": `Participant-aware modelling in ${c.episodes} selected episodes defines the current evidential boundary.`,
    },
    9: {
      "Participant 20026: symptoms emerge with a stable VF profile": `${c.toUti} operational Not_UTI-to-UTI transitions are all genomically linked`,
      "A concrete transition clarifies what same-strain, stable-profile evidence means.": `${c.toUtiAtThreshold} of ${c.toUti} focused direct pairs are at or below ${c.snpThreshold} SNPs.`,
      "T3": `${c.toUti} cases`,
      "culture-supported episode\nwithout UTI classification": "operational Not_UTI\nstarting endpoints",
      "UTI-1": `${c.casebookLinked} linked`,
      "symptom emergence on\nculture-supported bacteriuria": "operational UTI\nending endpoints",
      "42 DAYS": `${c.snpThreshold} SNP CUTOFF`,
      "5 SNPs": `${c.toUtiAtThreshold} of ${c.toUti}`,
      "strong same-strain evidence": "at or below the operational threshold",
      "+0 / -0": `${c.transitionsAtThreshold}/${c.transitions}`,
      "stable": "adjacent pairs <= threshold",
      "5": String(c.directPairs),
      "low": "all direct pairs",
      "Not_UTI -> UTI": `${c.toUti} -> ${c.casebookLinked}`,
      "changed": "cases -> linked",
      "Supports a host-state or unmeasured-regulation hypothesis; it does not prove mechanism.": "Close direct pairs motivate host-state or regulatory hypotheses; neither close nor distant pairs prove mechanism.",
    },
    10: {
      "Eleven clinical transitions include ten with linked WGS/VF evidence and one missing endpoint.": `${c.casebookCases} operational transitions are all genomically linked; ${c.casebookMissing} endpoints are missing.`,
    },
    11: {
      "Stable measured profiles can accompany changing clinical state": "Direct paired evidence distinguishes close from distant clinical transitions",
      "Four transitions show same-strain stable profiles; three are consistent with strain replacement.": `${c.toUtiAtThreshold}/${c.toUti} focused transitions and ${c.transitionsAtThreshold}/${c.transitions} adjacent pairs are <=${c.snpThreshold} SNPs.`,
    },
    12: {
      "No global VF association remains significant after FDR correction": "Population-level VF analysis remains exploratory",
      "With only 17 primary UTI VF-ready episodes, global signals remain exploratory.": `With ${c.uti} operational UTI episodes, participant-aware VF signals remain exploratory and non-causal.`,
    },
    13: {
      "Stable VF profiles occur across some symptom-emergence transitions.": `${c.toUtiAtThreshold} of ${c.toUti} focused direct pairs are <=${c.snpThreshold} SNPs.`,
      "Replacement is also present; no single mechanism fits all cases.": `${c.toUti - c.toUtiAtThreshold} focused pairs are above the threshold; no single mechanism is established.`,
      "Only 17 UTI episodes are VF/model-ready.": `Only ${c.uti} selected episodes meet the operational UTI phenotype.`,
    },
    14: {
      "Canonical outputs and rerun path": "Longcycler release outputs and rerun path",
      "Use these files as current sources of truth when navigating or reproducing the review.": `Claim registry -> ${c.directPairs} direct pairs -> ${c.transitions} adjacent pairs -> RQ01-RQ10.`,
      "results/summary/final_key_results_summary.md": "results/pipeline/longcycler_release_claim_registry.json",
      "results/final_figures/final_figure_captions.md": "results/clinical/analysis_cohort_longcycler.csv",
      "results/final_figures/final_figure_manifest.csv": "results/strain_compare/pairwise_metrics.csv",
      "results/mechanism/not_uti_to_uti_casebook.md": "results/mechanism/not_uti_to_uti_casebook.csv",
      "docs/pipeline_architecture.md": "results/research_questions/RUN_COMPLETE.txt",
      "01  Validate upstream outputs": "01  Validate selected cohort",
      "02  Confirm status_map.csv": "02  Confirm selected manifest",
      "03  Confirm vf_analysis_ready.csv": "03  Confirm direct pair evidence",
      "04  Refresh casebook": "04  Run RQ01-RQ10",
      "05  Render final figures": "05  Publish claim registry and decks",
      "35_final_figure_pack.R": "RUN_COMPLETE_ANALYSIS.sh",
      "Keep primary UTI_Status outputs separate from superseded legacy labels.": "Use the versioned operational phenotype and selected Longcycler cohort throughout.",
    },
    15: {
      "Current UTI classification audit": "Operational UTI phenotype audit",
      "Current catheter-aware symptom and culture-support rule; diagnostic context only.": `Versioned catheter-aware symptom and culture-support rule applied to ${c.episodes} selected episodes.`,
    },
    19: {
      "Near-miss episodes remain sensitivity-only and do not enter the primary UTI group.": `${c.nearMiss} near-miss rows remain sensitivity-only and do not enter the operational UTI group.`,
    },
    20: {
      "Instability is expected with 17 VF-ready UTI episodes and is not association testing.": `Instability is expected with ${c.uti} operational UTI episodes and is not confirmatory association testing.`,
    },
    22: {
      "-  583 clinical episodes: 18 UTI / 565 Not_UTI": `-  ${c.episodes} selected episodes: ${c.uti} UTI / ${c.notUti} Not_UTI`,
      "-  556 VF-ready episodes: 17 UTI / 539 Not_UTI": `-  ${c.directPairs} direct pairs; ${c.transitions} adjacent from ${c.transitionResidents} residents`,
      "-  11 clinical transitions; 10 WGS/VF-linked": `-  ${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing} focused cases / linked / missing; ${c.toUtiAtThreshold} <=${c.snpThreshold} SNPs`,
      "-  No global VF association significant after FDR": `-  ${c.nearMiss} near-miss rows are sensitivity-only; RQ01-RQ10 complete`,
      "-  ASB-versus-UTI as the current comparison": "-  a non-operational phenotype as the current comparison",
      "-  significant lpf findings as current evidence": "-  named VF or mutation signals as causal evidence",
      "results/longitudinal/swimmer_plot.png  |  full longitudinal context": "plots/publication/Fig1_Swimmer_Plot.png  |  contextual longitudinal lookup",
    },
  };
}

function deckTextMap(deckKey, counts) {
  if (deckKey === "onboarding") return onboardingTextMap(counts);
  if (deckKey === "scientific") return scientificTextMap(counts);
  if (deckKey === "vf_focused") return vfTextMap(counts);
  throw new Error(`Unknown deck key: ${deckKey}`);
}

function resolvePlotPath(registry, suffix) {
  const normalizedSuffix = suffix.replaceAll("\\", "/");
  const matches = (registry.plot_files || []).filter((candidate) =>
    String(candidate).replaceAll("\\", "/").endsWith(normalizedSuffix),
  );
  if (matches.length !== 1) {
    throw new Error(`Expected exactly one current plot for ${suffix}; found ${matches.length}`);
  }
  const resolved = path.resolve(matches[0]);
  if (resolved.toLowerCase().includes(RETIRED_INPUT_TOKEN)) {
    throw new Error(`Retired input token in plot path: ${resolved}`);
  }
  return resolved;
}

async function makeFixturePlot(artifact, outputPath, label, c, accentIndex) {
  const { Presentation } = artifact;
  const p = Presentation.create({ slideSize: { width: 1280, height: 720 } });
  const slide = p.slides.add();
  slide.background.fill = "white";
  const palette = ["#2C7FB8", "#2E8B62", "#D95F0E", "#6A51A3", "#4F6D7A"];
  const accent = palette[accentIndex % palette.length];
  const title = slide.shapes.add({
    geometry: "textbox",
    name: "fixture-title",
    position: { left: 62, top: 42, width: 1156, height: 72 },
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  title.text = label;
  title.text.style = { fontSize: 30, bold: true, color: "#172B36", typeface: "Arial" };
  const subtitle = slide.shapes.add({
    geometry: "textbox",
    name: "fixture-subtitle",
    position: { left: 62, top: 110, width: 1156, height: 48 },
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  subtitle.text = "Development fixture — final build resolves the current plot from the release registry";
  subtitle.text.style = { fontSize: 18, color: "#53636A", typeface: "Arial" };

  const metrics = [
    ["Selected episodes", c.episodes, c.episodes],
    ["Direct pairs", c.directPairs, c.directPairs],
    ["Adjacent pairs", c.transitions, c.directPairs],
    ["Not_UTI to UTI", c.toUti, c.directPairs],
  ];
  metrics.forEach(([name, value, max], index) => {
    const y = 202 + index * 98;
    const labelShape = slide.shapes.add({
      geometry: "textbox",
      name: `fixture-label-${index + 1}`,
      position: { left: 74, top: y, width: 260, height: 38 },
      fill: "none",
      line: { style: "solid", fill: "none", width: 0 },
    });
    labelShape.text = `${name}  ${value}`;
    labelShape.text.style = { fontSize: 18, bold: true, color: "#172B36", typeface: "Arial" };
    slide.shapes.add({
      geometry: "rect",
      name: `fixture-track-${index + 1}`,
      position: { left: 344, top: y + 4, width: 800, height: 28 },
      fill: "#E7ECEF",
      line: { style: "solid", fill: "none", width: 0 },
    });
    slide.shapes.add({
      geometry: "rect",
      name: `fixture-bar-${index + 1}`,
      position: { left: 344, top: y + 4, width: Math.max(18, 800 * Number(value) / Number(max)), height: 28 },
      fill: accent,
      line: { style: "solid", fill: "none", width: 0 },
    });
  });
  const footer = slide.shapes.add({
    geometry: "textbox",
    name: "fixture-footer",
    position: { left: 62, top: 624, width: 1156, height: 48 },
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  footer.text = `${c.episodes}/${c.residents}/${c.uti}/${c.notUti} | ${c.transitions}/${c.transitionResidents} | casebook ${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing}`;
  footer.text.style = { fontSize: 17, color: "#53636A", typeface: "Arial" };
  await writeBlob(outputPath, await p.export({ slide, format: "png", scale: 1.25 }));
  return outputPath;
}

async function replaceImage(image, bytes, alt) {
  const frame = image.frame;
  const crop = image.crop;
  const geometry = image.geometry;
  const borderRadius = image.borderRadius;
  const rotation = image.rotation;
  const flipHorizontal = image.flipHorizontal;
  const flipVertical = image.flipVertical;
  const lockAspectRatio = image.lockAspectRatio;
  image.replace({ blob: fileBytes(bytes), contentType: "image/png", alt, fit: "contain" });
  image.frame = frame;
  image.crop = crop;
  image.geometry = geometry;
  image.borderRadius = borderRadius;
  image.rotation = rotation;
  image.flipHorizontal = flipHorizontal;
  image.flipVertical = flipVertical;
  image.lockAspectRatio = lockAspectRatio;
}

function addWorkflowBox(slide, name, position, fill, line, title, body) {
  const box = slide.shapes.add({
    geometry: "roundRect",
    name,
    position,
    fill,
    line: { style: "solid", fill: line, width: 1.4 },
    borderRadius: 8,
  });
  box.text = `${title}\n${body}`;
  box.text.style = { fontSize: 16, color: "#172B36", typeface: "Arial", bold: false };
  return box;
}

function addWorkflowArrow(slide, name, position) {
  return slide.shapes.add({
    geometry: "rightArrow",
    name,
    position,
    fill: "#7B8B95",
    line: { style: "solid", fill: "none", width: 0 },
  });
}

function addScientificWorkflow(slide, slideNumber, c) {
  if (slideNumber === 6) {
    [
      { left: 275, top: 278, width: 15, height: 24 },
      { left: 472, top: 278, width: 15, height: 24 },
      { left: 669, top: 278, width: 15, height: 24 },
    ].forEach((position, index) => {
      addWorkflowArrow(slide, `lc-native-connector-${index + 1}`, position);
    });
    [
      addWorkflowBox(slide, "lc-native-cohort", { left: 95, top: 215, width: 178, height: 150 }, "#E8F2F8", "#2C7FB8", "SELECTED COHORT", `${c.episodes} episodes\n${c.residents} residents\n${c.uti}/${c.notUti} UTI/Not_UTI`),
      addWorkflowBox(slide, "lc-native-genomes", { left: 292, top: 215, width: 178, height: 150 }, "#E8F4EE", "#2E8B62", "SELECTED GENOMES", `${c.episodes} QC-passing\nLongcycler assemblies\none per episode`),
      addWorkflowBox(slide, "lc-native-pairs", { left: 489, top: 215, width: 178, height: 150 }, "#EEF1F4", "#6B7C8F", "DIRECT EVIDENCE", `${c.directPairs} resident pairs\n${c.transitions} adjacent\n${c.transitionsAtThreshold} <=${c.snpThreshold} SNPs`),
      addWorkflowBox(slide, "lc-native-cases", { left: 686, top: 215, width: 178, height: 150 }, "#FBEDE5", "#D95F0E", "FOCUSED CASEBOOK", `${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing}\ncases/linked/missing\n${c.toUtiAtThreshold}/${c.toUti} <=${c.snpThreshold} SNPs`),
    ];
    const banner = slide.shapes.add({
      geometry: "roundRect",
      name: "lc-native-release-banner",
      position: { left: 95, top: 400, width: 769, height: 80 },
      fill: "#F4F6F7",
      line: { style: "solid", fill: "#D3DCE1", width: 1 },
      borderRadius: 7,
    });
    banner.text = `Release layer: operational phenotype -> direct pairs -> adjacent transitions -> RQ01-RQ10`;
    banner.text.style = { fontSize: 17, bold: true, color: "#172B36", typeface: "Arial" };
    return;
  }
  if (slideNumber === 14) {
    [
      { left: 334, top: 310, width: 27, height: 24 },
      { left: 600, top: 310, width: 26, height: 24 },
    ].forEach((position, index) => {
      addWorkflowArrow(slide, `lc-native-connector-${index + 1}`, position);
    });
    [
      addWorkflowBox(slide, "lc-native-registry", { left: 105, top: 250, width: 220, height: 145 }, "#E8F2F8", "#2C7FB8", "1 | CLAIM REGISTRY", "Validated counts\nplot inventory\nsource hashes"),
      addWorkflowBox(slide, "lc-native-evidence", { left: 370, top: 250, width: 220, height: 145 }, "#E8F4EE", "#2E8B62", "2 | EVIDENCE TABLES", `${c.directPairs} direct pairs\n${c.transitions} adjacent pairs\n${c.casebookCases} linked cases`),
      addWorkflowBox(slide, "lc-native-release", { left: 635, top: 250, width: 220, height: 145 }, "#FBEDE5", "#D95F0E", "3 | RELEASE OUTPUTS", "RQ01-RQ10\ncurrent final plots\nregenerated decks"),
    ];
    const runner = slide.shapes.add({
      geometry: "roundRect",
      name: "lc-native-runner",
      position: { left: 105, top: 435, width: 750, height: 70 },
      fill: "#F4F6F7",
      line: { style: "solid", fill: "#D3DCE1", width: 1 },
      borderRadius: 7,
    });
    runner.text = "Single entry point: RUN_COMPLETE_ANALYSIS.sh";
    runner.text.style = { fontSize: 18, bold: true, color: "#172B36", typeface: "Arial" };
  }
}

function nativeWorkflowNames(slideNumber) {
  if (slideNumber === 6) {
    return [
      "lc-native-connector-1", "lc-native-connector-2", "lc-native-connector-3",
      "lc-native-cohort", "lc-native-genomes", "lc-native-pairs", "lc-native-cases",
      "lc-native-release-banner",
    ];
  }
  if (slideNumber === 14) {
    return [
      "lc-native-connector-1", "lc-native-connector-2",
      "lc-native-registry", "lc-native-evidence", "lc-native-release", "lc-native-runner",
    ];
  }
  return [];
}

function nativeWorkflowText(slideNumber, c) {
  if (slideNumber === 6) {
    return new Map([
      ["lc-native-cohort", `SELECTED COHORT\n${c.episodes} episodes\n${c.residents} residents\n${c.uti}/${c.notUti} UTI/Not_UTI`],
      ["lc-native-genomes", `SELECTED GENOMES\n${c.episodes} QC-passing\nLongcycler assemblies\none per episode`],
      ["lc-native-pairs", `DIRECT EVIDENCE\n${c.directPairs} resident pairs\n${c.transitions} adjacent\n${c.transitionsAtThreshold} <=${c.snpThreshold} SNPs`],
      ["lc-native-cases", `FOCUSED CASEBOOK\n${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing}\ncases/linked/missing\n${c.toUtiAtThreshold}/${c.toUti} <=${c.snpThreshold} SNPs`],
      ["lc-native-release-banner", "Release layer: operational phenotype -> direct pairs -> adjacent transitions -> RQ01-RQ10"],
    ]);
  }
  if (slideNumber === 14) {
    return new Map([
      ["lc-native-registry", "1 | CLAIM REGISTRY\nValidated counts\nplot inventory\nsource hashes"],
      ["lc-native-evidence", `2 | EVIDENCE TABLES\n${c.directPairs} direct pairs\n${c.transitions} adjacent pairs\n${c.casebookCases} linked cases`],
      ["lc-native-release", "3 | RELEASE OUTPUTS\nRQ01-RQ10\ncurrent final plots\nregenerated decks"],
      ["lc-native-runner", "Single entry point: RUN_COMPLETE_ANALYSIS.sh"],
    ]);
  }
  return new Map();
}

function refreshOrAddScientificWorkflow({ presentation, slide, slideNumber, counts, importedLayout }) {
  const expectedNames = nativeWorkflowNames(slideNumber);
  const existing = (importedLayout.elements || []).filter(
    (element) => String(element.name || "").startsWith("lc-native-"),
  );
  if (!existing.length) {
    addScientificWorkflow(slide, slideNumber, counts);
    return "added";
  }
  const existingNames = existing.map((element) => element.name);
  const duplicateNames = existingNames.filter((name, index) => existingNames.indexOf(name) !== index);
  const missingNames = expectedNames.filter((name) => !existingNames.includes(name));
  const unexpectedNames = existingNames.filter((name) => !expectedNames.includes(name));
  if (duplicateNames.length || missingNames.length || unexpectedNames.length || existing.length !== expectedNames.length) {
    throw new Error(
      `scientific slide ${slideNumber}: non-idempotent native workflow state `
      + JSON.stringify({ duplicateNames, missingNames, unexpectedNames, observed: existingNames }),
    );
  }
  for (const [name, nextText] of nativeWorkflowText(slideNumber, counts)) {
    const element = existing.find((candidate) => candidate.name === name);
    const target = presentation.resolve(element.aid);
    if (!target?.text) throw new Error(`scientific slide ${slideNumber}: could not refresh ${name}`);
    target.text.set(nextText);
  }
  return "refreshed";
}

async function prepareStarter(projectRoot, workspace, deck, templatePptx, registryPath) {
  await fs.rm(workspace, { recursive: true, force: true });
  await fs.mkdir(workspace, { recursive: true });
  await runTool(process.execPath, [
    path.join(SKILL_DIR, "template_following_scripts", "inspect_template_deck.mjs"),
    "--workspace", workspace,
    "--pptx", templatePptx,
    "--scale", "1.5",
  ]);
  await runTool(process.execPath, [
    path.join(projectRoot, "scripts", "generated", "create_core_review_template_map.mjs"),
    "--workspace", workspace,
    "--deck-kind", deck.key,
    "--registry", registryPath,
  ]);
  await runTool(process.execPath, [
    path.join(SKILL_DIR, "template_following_scripts", "validate_template_plan.mjs"),
    "--workspace", workspace,
    "--map", path.join(workspace, "template-frame-map.json"),
  ], { echo: true });
  await runTool(process.execPath, [
    path.join(SKILL_DIR, "template_following_scripts", "prepare_template_starter_deck.mjs"),
    "--workspace", workspace,
    "--pptx", templatePptx,
    "--map", path.join(workspace, "template-frame-map.json"),
    "--out", path.join(workspace, "template-starter.pptx"),
    "--preview-dir", path.join(workspace, "template-starter-preview"),
    "--layout-dir", path.join(workspace, "template-starter-layout"),
    "--scale", "1.5",
  ], { env: { PYTHON }, echo: true });
  await runTool(PYTHON, [
    path.join(SKILL_DIR, "container_tools", "create_montage.py"),
    "--input_dir", path.join(workspace, "template-starter-preview"),
    "--output_file", path.join(workspace, "template-starter-contact-sheet.png"),
    "--num_col", "4",
    "--label_mode", "filename",
    "--fail_on_image_error",
  ]);
}

function assertVisibleText(deck, text, c) {
  const lower = text.toLowerCase();
  const stalePatterns = [
    ["583", /\b583\b/i],
    ["585", /\b585\b/i],
    ["565", /\b565\b/i],
    ["556", /\b556\b/i],
    ["539", /\b539\b/i],
    ["394", /\b394\b/i],
    ["162", /\b162\b/i],
    ["116", /\b116\b/i],
    ["7/9", /\b7\s*(?:\/|of)\s*9\b/i],
    ["1,291", /\b1,291\b/i],
    ["11 clinical", /\b11\s+clinical\b/i],
    ["10 WGS/VF", /\b10\s+wgs\/vf\b/i],
    ["17 UTI", /\b17\s+(?:operational\s+)?uti\b/i],
    ["RQ11", /\brq11\b/i],
  ];
  for (const [label, pattern] of stalePatterns) {
    if (pattern.test(text)) throw new Error(`${deck.key}: stale visible text remains: ${label}`);
  }
  if (lower.includes(RETIRED_INPUT_TOKEN)) throw new Error(`${deck.key}: retired input token remains in visible text`);
  for (const required of [
    String(c.episodes), String(c.residents), String(c.uti), String(c.notUti),
    String(c.directPairs), String(c.transitions), String(c.transitionResidents),
    String(c.toUti), String(c.toUtiAtThreshold), String(c.nearMiss),
  ]) {
    if (!new RegExp(`\\b${required}\\b`).test(lower)) {
      throw new Error(`${deck.key}: required release anchor missing from visible text: ${required}`);
    }
  }
  if (!lower.includes("rq01-rq10")) throw new Error(`${deck.key}: required release anchor missing from visible text: RQ01-RQ10`);
  if (!lower.includes("operational")) throw new Error(`${deck.key}: operational phenotype label missing`);
}

function plannedImageSlides(deck) {
  return new Set([
    ...Object.keys(deck.imagePlan || {}).map(Number),
    ...Array.from(deck.nativeSlides || []),
  ]);
}

function validateMappedImageAuthorizations(deck, frameMap) {
  const expected = plannedImageSlides(deck);
  const bySlide = new Map();
  for (const entry of frameMap.outputSlides || []) {
    if (bySlide.has(entry.outputSlide)) {
      throw new Error(`${deck.key}: duplicate frame-map entry for slide ${entry.outputSlide}`);
    }
    bySlide.set(entry.outputSlide, entry);
    const replacements = (entry.editTargets || []).filter((target) => target.action === "replace");
    if (!expected.has(entry.outputSlide)) {
      if (replacements.length) {
        throw new Error(`${deck.key} slide ${entry.outputSlide}: inherited image replacement is not in the deck image plan`);
      }
      continue;
    }
    if (replacements.length !== 1) {
      throw new Error(`${deck.key} slide ${entry.outputSlide}: expected one mapped image authorization; found ${replacements.length}`);
    }
    const target = replacements[0];
    if (typeof target.sourceElementId !== "string" || !target.sourceElementId || target.sourceElementIds !== undefined) {
      throw new Error(`${deck.key} slide ${entry.outputSlide}: mapped image authorization must name one exact sourceElementId`);
    }
  }
  for (const slideNumber of expected) {
    if (!bySlide.has(slideNumber)) throw new Error(`${deck.key} slide ${slideNumber}: missing frame-map entry`);
  }
  return bySlide;
}

function bboxesEquivalent(left, right, tolerance = 0.75) {
  return Array.isArray(left)
    && Array.isArray(right)
    && left.length === 4
    && right.length === 4
    && left.every((value, index) => Math.abs(Number(value) - Number(right[index])) <= tolerance);
}

async function resolveMappedStarterImageElement({ workspace, deck, slideNumber, starterLayout, frameMapBySlide }) {
  const frameEntry = frameMapBySlide.get(slideNumber);
  const replacement = (frameEntry?.editTargets || []).find((target) => target.action === "replace");
  if (!replacement?.sourceElementId) {
    throw new Error(`${deck.key} slide ${slideNumber}: exact mapped source image ID is unavailable`);
  }
  const padded = String(slideNumber).padStart(2, "0");
  const sourceLayout = JSON.parse(
    await fs.readFile(
      path.join(workspace, "template-inspect", "layouts", `source-slide-${padded}.layout.json`),
      "utf8",
    ),
  );
  const sourceMatches = (sourceLayout.elements || []).filter(
    (element) => element.kind === "image" && element.aid === replacement.sourceElementId,
  );
  if (sourceMatches.length !== 1 || !sourceMatches[0].bbox) {
    throw new Error(`${deck.key} slide ${slideNumber}: mapped source image ${replacement.sourceElementId} did not resolve uniquely`);
  }
  const starterMatches = (starterLayout.elements || []).filter(
    (element) => element.kind === "image" && bboxesEquivalent(element.bbox, sourceMatches[0].bbox),
  );
  if (starterMatches.length !== 1 || !starterMatches[0].aid) {
    throw new Error(
      `${deck.key} slide ${slideNumber}: mapped source image ${replacement.sourceElementId} matched ${starterMatches.length} starter images`,
    );
  }
  return {
    sourceAid: replacement.sourceElementId,
    starterAid: starterMatches[0].aid,
    starterElement: starterMatches[0],
  };
}

function bboxInside(inner, outer, tolerance = 1.5) {
  const [il, it, iw, ih] = inner;
  const [ol, ot, ow, oh] = outer;
  return il >= ol - tolerance && it >= ot - tolerance && il + iw <= ol + ow + tolerance && it + ih <= ot + oh + tolerance;
}

function intersectionArea(a, b) {
  const left = Math.max(a[0], b[0]);
  const top = Math.max(a[1], b[1]);
  const right = Math.min(a[0] + a[2], b[0] + b[2]);
  const bottom = Math.min(a[1] + a[3], b[1] + b[3]);
  return Math.max(0, right - left) * Math.max(0, bottom - top);
}

async function auditNativeLayouts(deck, workspace, frameMapBySlide) {
  const results = [];
  for (const slideNumber of deck.nativeSlides || []) {
    const padded = String(slideNumber).padStart(2, "0");
    const starter = JSON.parse(await fs.readFile(path.join(workspace, "template-starter-layout", `starter-slide-${padded}.layout.json`), "utf8"));
    const final = JSON.parse(await fs.readFile(path.join(workspace, "layout", "final", `slide-${padded}.layout.json`), "utf8"));
    const mappedImage = await resolveMappedStarterImageElement({
      workspace,
      deck,
      slideNumber,
      starterLayout: starter,
      frameMapBySlide,
    });
    const frame = mappedImage.starterElement.bbox;
    if (!frame) throw new Error(`${deck.key} slide ${slideNumber}: missing inherited image frame for native-layout audit`);
    const added = final.elements.filter((element) => String(element.name || "").startsWith("lc-native-") && element.bbox);
    if (!added.length) throw new Error(`${deck.key} slide ${slideNumber}: native rebuild elements not found`);
    const duplicateNames = added.map((element) => element.name)
      .filter((name, index, names) => names.indexOf(name) !== index);
    if (duplicateNames.length) {
      throw new Error(`${deck.key} slide ${slideNumber}: duplicate native-shape names: ${JSON.stringify(duplicateNames)}`);
    }
    for (const element of added) {
      if (!bboxInside(element.bbox, frame)) {
        throw new Error(`${deck.key} slide ${slideNumber}: ${element.name} escapes the inherited image frame`);
      }
    }
    const cards = added.filter((element) => !String(element.name).includes("connector"));
    const overlaps = [];
    for (let i = 0; i < cards.length; i += 1) {
      for (let j = i + 1; j < cards.length; j += 1) {
        const area = intersectionArea(cards[i].bbox, cards[j].bbox);
        if (area > 4) overlaps.push([cards[i].name, cards[j].name, area]);
      }
    }
    if (overlaps.length) throw new Error(`${deck.key} slide ${slideNumber}: unintended native-shape overlaps: ${JSON.stringify(overlaps)}`);
    results.push({ slide: slideNumber, inherited_frame: frame, native_elements: added.map((element) => element.name), status: "PASS" });
  }
  await writeJsonAtomic(path.join(workspace, "native_overlap_audit.json"), { deck: deck.key, results });
}

async function buildDeck({ artifact, projectRoot, baseWorkspace, outputRoot, deck, registry, registryPath, counts, date, devFixture }) {
  const workspace = path.join(baseWorkspace, deck.key);
  const sourcePptx = path.join(projectRoot, deck.source);
  const outputPptx = outputRoot ? path.join(outputRoot, deck.output) : path.join(projectRoot, deck.output);
  await prepareStarter(projectRoot, workspace, deck, sourcePptx, registryPath);
  await fs.writeFile(
    path.join(workspace, "authoring-source.mjs"),
    `${buildDeck.toString()}\n`,
    "utf8",
  );
  const frameMap = JSON.parse(await fs.readFile(path.join(workspace, "template-frame-map.json"), "utf8"));
  const frameMapBySlide = validateMappedImageAuthorizations(deck, frameMap);

  const { FileBlob, PresentationFile } = artifact;
  const presentation = await PresentationFile.importPptx(await FileBlob.load(path.join(workspace, "template-starter.pptx")));
  const slides = slidesFromPresentation(presentation);
  if (slides.length !== deck.slides) throw new Error(`${deck.key}: expected ${deck.slides} slides; imported ${slides.length}`);
  const textMap = deckTextMap(deck.key, counts);
  const previewDir = path.join(workspace, "final-preview");
  const layoutDir = path.join(workspace, "layout", "final");
  const assetDir = path.join(workspace, "assets");
  await fs.mkdir(previewDir, { recursive: true });
  await fs.mkdir(layoutDir, { recursive: true });
  await fs.mkdir(assetDir, { recursive: true });
  const provenance = [];
  const visibleText = [];

  for (let slideNumber = 1; slideNumber <= slides.length; slideNumber += 1) {
    const slide = slides[slideNumber - 1];
    const padded = String(slideNumber).padStart(2, "0");
    // Resolve editable objects from the current in-memory import. Artifact IDs
    // captured while preparing template-starter.pptx may be regenerated by a
    // subsequent import even though the inherited objects and bounds persist.
    const importedLayout = JSON.parse(await (await slide.export({ format: "layout" })).text());
    const textElements = importedLayout.elements.filter((element) => typeof element.text === "string");
    const sourceOccurrences = new Map();
    for (const element of textElements) {
      const target = presentation.resolve(element.aid);
      if (!target?.text) throw new Error(`${deck.key} slide ${slideNumber}: unresolved inherited text anchor ${element.aid}`);
      const occurrenceIndex = sourceOccurrences.get(element.text) || 0;
      sourceOccurrences.set(element.text, occurrenceIndex + 1);
      const next = lookup(textMap, slideNumber, element.text, date, element.aid, occurrenceIndex);
      target.text.set(next);
      const finalText = applyInheritedTextVisualRepair({
        deckKey: deck.key,
        slideNumber,
        target,
        element,
        next,
      });
      visibleText.push(finalText);
    }

    let mappedImageTarget = null;
    if (deck.nativeSlides?.has(slideNumber) || deck.imagePlan[slideNumber]) {
      const mappedImage = await resolveMappedStarterImageElement({
        workspace,
        deck,
        slideNumber,
        starterLayout: importedLayout,
        frameMapBySlide,
      });
      mappedImageTarget = presentation.resolve(mappedImage.starterAid);
      if (!mappedImageTarget || typeof mappedImageTarget.replace !== "function") {
        throw new Error(`${deck.key} slide ${slideNumber}: mapped starter image ${mappedImage.starterAid} is not replaceable`);
      }
    }

    if (deck.nativeSlides?.has(slideNumber)) {
      await replaceImage(mappedImageTarget, TRANSPARENT_PNG, "Transparent replacement for editable native workflow rebuild");
      refreshOrAddScientificWorkflow({ presentation, slide, slideNumber, counts, importedLayout });
      provenance.push({ deck: deck.key, slide: slideNumber, asset: "editable native PowerPoint workflow", sha256: "native-shapes", source: registryPath });
    } else if (deck.imagePlan[slideNumber]) {
      let assetPath;
      if (devFixture) {
        assetPath = path.join(assetDir, `fixture-${padded}.png`);
        await makeFixturePlot(artifact, assetPath, path.basename(deck.imagePlan[slideNumber], ".png").replaceAll("_", " "), counts, slideNumber);
      } else {
        assetPath = resolvePlotPath(registry, deck.imagePlan[slideNumber]);
      }
      const bytes = await fs.readFile(assetPath);
      if (!bytes.length) throw new Error(`${deck.key} slide ${slideNumber}: empty replacement asset ${assetPath}`);
      await replaceImage(mappedImageTarget, bytes, `Current final-run plot: ${path.basename(assetPath)}`);
      provenance.push({ deck: deck.key, slide: slideNumber, asset: path.basename(assetPath), sha256: hashBuffer(bytes), source: assetPath });
    }

    slide.speakerNotes.textFrame.setText([
      `Release source: ${registryPath}`,
      "Scope: selected QC-passing Longcycler assemblies only.",
      "Clinical phenotype: versioned operational UTI phenotype.",
      `Analytical anchors: ${counts.episodes} episodes; ${counts.directPairs} direct pairs; ${counts.transitions} adjacent pairs; RQ01-RQ10.`,
    ].join("\n"));
    slide.speakerNotes.setVisible(true);

    await writeBlob(path.join(previewDir, `slide-${padded}.png`), await presentation.export({ slide, format: "png", scale: 2 }));
    await writeBlob(path.join(layoutDir, `slide-${padded}.layout.json`), await slide.export({ format: "layout" }));
  }

  assertVisibleText(deck, visibleText.join("\n"), counts);
  await writeCsvAtomic(path.join(workspace, "image-provenance.csv"), provenance, ["deck", "slide", "asset", "sha256", "source"]);
  await writeBlob(path.join(workspace, "final-montage.webp"), await presentation.export({ format: "webp", montage: true, scale: 1 }));
  const inspect = await presentation.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 500000 });
  await fs.writeFile(path.join(workspace, "final-inspect.ndjson"), inspect.ndjson || "", "utf8");

  await fs.mkdir(path.dirname(outputPptx), { recursive: true });
  const tempPptx = path.join(workspace, `${deck.key}-final.pptx`);
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(tempPptx);
  await assertNoEmptyPlaceholders(tempPptx);
  await assertPptxHasNoForbiddenContent(tempPptx);
  await auditNativeLayouts(deck, workspace, frameMapBySlide);
  await runTool(process.execPath, [
    path.join(SKILL_DIR, "template_following_scripts", "check_template_fidelity.mjs"),
    "--workspace", workspace,
    "--starter-pptx", path.join(workspace, "template-starter.pptx"),
    "--final-pptx", tempPptx,
    "--map", path.join(workspace, "template-frame-map.json"),
    "--starter-layout-dir", path.join(workspace, "template-starter-layout"),
    "--final-layout-dir", layoutDir,
    "--edit-dir", workspace,
  ], { echo: true });
  await auditOverflowAgainstStarter({
    deck,
    workspace,
    starterPptx: path.join(workspace, "template-starter.pptx"),
    finalPptx: tempPptx,
  });

  const targetTemp = path.join(path.dirname(outputPptx), `.${path.basename(outputPptx)}.${process.pid}.tmp.pptx`);
  await fs.copyFile(tempPptx, targetTemp);
  await fs.rename(targetTemp, outputPptx);
  return { deck: deck.key, outputPptx, workspace, slides: deck.slides, imageReplacements: provenance.length };
}

async function main() {
  if (await bootstrapToScratch({
    scriptUrl: import.meta.url,
    taskSlug: "ruti-core-longcycler-review-decks",
    argv: process.argv.slice(2),
  })) return;

  const args = parseArgs(process.argv.slice(2));
  const projectRoot = path.resolve(args["project-root"] || DEFAULT_PROJECT_ROOT);
  const baseWorkspace = path.resolve(args.workspace);
  const devFixture = Boolean(args["dev-fixture"]);
  const outputRoot = args["output-root"]
    ? path.resolve(args["output-root"])
    : devFixture
      ? path.join(baseWorkspace, "development-output")
      : null;
  const { registry, registryPath } = await loadRegistry({
    projectRoot,
    registryPath: args.registry,
    devFixture,
  });
  const counts = releaseCounts(registry);
  const date = releaseDate(registry);
  const artifact = await import("@oai/artifact-tool");
  const results = [];
  for (const deck of DECKS) {
    results.push(await buildDeck({
      artifact,
      projectRoot,
      baseWorkspace,
      outputRoot,
      deck,
      registry,
      registryPath,
      counts,
      date,
      devFixture,
    }));
  }
  await writeJsonAtomic(path.join(baseWorkspace, "core-review-deck-build-summary.json"), {
    mode: devFixture ? "development fixture" : "final registry",
    registry: registryPath,
    results,
  });
  console.log(JSON.stringify(results, null, 2));
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

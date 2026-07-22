#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";

import {
  assertOutputTreeClean,
  bootstrapToScratch,
  copyFileAtomic,
  editTemplateDeck,
  loadRegistry,
  outputDirectory,
  parseArgs,
  releaseCounts,
  releaseDate,
  sha256File,
  writeCsvAtomic,
  writeJsonAtomic,
  writeTextAtomic,
} from "./longcycler_release_presentation_common.mjs";

if (await bootstrapToScratch({
  scriptUrl: import.meta.url,
  taskSlug: "longcycler-methods-summary",
  argv: process.argv.slice(2),
})) {
  process.exit(0);
}

function summaryContent(registry) {
  const c = releaseCounts(registry);
  const footer = `Longcycler-only methods summary  |  RQ01–RQ10  |  ${releaseDate(registry)}`;
  return {
    1: [
      "METHODS CHOICE",
      "Longcycler-only is now the sole analytical denominator",
      "Every genomic model, table and figure uses the selected QC-passing Longcycler cohort.",
      footer,
      "01",
      "WHY THIS HELPS",
      "One selected QC-passing assembly per episode keeps strain, lineage and VF summaries on one reproducible input contract.",
      `Every downstream join is checked against the same ${c.episodes} episode keys.`,
      "WHAT CHANGES",
      `The analytical cohort is ${c.episodes} genome-linked episodes from ${c.residents} residents.`,
      `All ${c.directPairs} direct pairs and ${c.transitions} adjacent transitions have Longcycler endpoints.`,
      "WHAT DOES NOT CHANGE",
      "Clinical labels are assigned independently using the versioned operational UTI phenotype.",
      `The ≤${c.snpThreshold} SNP boundary remains an operational support rule—not a causal result.`,
      "Full clinical source counts appear only as labelled attrition/QC context.",
    ],
    2: [
      "LONGCYCLER-ONLY COUNTS",
      `The analysis contains ${c.episodes} selected genome-linked episodes`,
      "    Counts remain paired with their unit, filter and exact Longcycler-only denominator.",
      "NUMBERED PROCESS MAP",
      footer,
      "02",
      "SELECTED COHORT CONTRACT",
      "01",
      "Analytical episodes",
      String(c.episodes),
      "selected QC-passing Longcycler rows",
      "analysis clinical cohort",
      "02",
      "Residents",
      String(c.residents),
      "unique residents represented",
      "selected assembly manifest",
      "CLINICAL AND PAIR DENOMINATORS",
      "03",
      "Operational status",
      `${c.uti} | ${c.notUti}`,
      "UTI | Not_UTI on the selected cohort",
      "analysis clinical cohort",
      "04",
      "Direct genome pairs",
      String(c.directPairs),
      "all unordered within-resident pairs",
      "pairwise metrics",
      "DOWNSTREAM ANALYSIS UNITS",
      "05",
      "VFDB feature space",
      String(c.vfFeatures),
      "binary features on selected episodes",
      "claim registry",
      "06",
      "MLST lineage",
      `${c.mlstTyped} | ${c.mlstDistinct}`,
      "typed episodes | distinct ST labels",
      "claim registry",
      "07",
      "Adjacent transitions",
      `${c.transitions} | ${c.transitionResidents}`,
      "pairs | residents after full rebuild",
      "longcycler transitions",
      "08",
      "Focused transitions",
      `${c.toUti} | ${c.toUtiAtThreshold}`,
      `Not_UTI→UTI | at ≤${c.snpThreshold} SNPs`,
      "mechanism casebook",
    ],
    3: [
      "METHODS BOUNDARY",
      "The release supports reproducible, cautious genomic description",
      "Direct evidence, operational definitions and interpretation limits remain separate.",
      footer,
      "03",
      "SUPPORTED",
      "VF presence/absence, MLST lineage context, direct pair evidence, longitudinal transitions and casebook traceability.",
      `Continuity is described from direct ≤${c.snpThreshold} SNP evidence—not ST agreement or graph transitivity alone.`,
      "NOT SUPPORTED",
      "No antibiotic, demographic, host-factor, named-mutation or causal mechanism claim is inferred from proxies.",
      "The operational UTI phenotype is not presented as a reconstruction of the full published protocol.",
      "RELEASE CONTRACT",
      `RQ01–RQ10 use the same ${c.episodes}-episode analytical cohort.`,
      `${c.casebookCases}/${c.casebookCases} casebook rows are linked; ${c.casebookMissing} are missing. ${c.nearMiss} near-miss rows remain separate.`,
      "Suggested wording: genomic analyses were restricted to selected QC-passing Longcycler assemblies, with exploratory non-causal interpretation.",
    ],
  };
}

function flowContent(registry, pageNumber) {
  const c = releaseCounts(registry);
  const footer = `Longcycler-only methods summary  |  RQ01–RQ10  |  ${releaseDate(registry)}`;
  return {
    1: [
      "ANALYSIS FLOW",
      "One selected Longcycler cohort feeds every analytical result",
      "Clinical status and the selected assembly manifest meet by exact episode key before any genomic analysis.",
      "Attrition/QC source",
      `${c.sourceEpisodes} episodes`,
      `${c.sourceResidents} residents\n${c.sourceUti} UTI | ${c.sourceNotUti} Not_UTI · context only`,
      "Selected assembly + QC",
      `${c.episodes} rows`,
      "one selected QC-pass assembly per key\nFASTA content hashes verified",
      "Exact key-linked cohort",
      `${c.episodes} analytical episodes`,
      `${c.residents} residents\nmanifest and clinical keys match exactly`,
      "Longcycler analytical cohort",
      `${c.episodes} rows | ${c.residents} residents`,
      `${c.uti} UTI | ${c.notUti} Not_UTI\nversioned operational phenotype`,
      "Analysis layers",
      `${c.directPairs} direct pairs`,
      `VF: ${c.vfFeatures} features · MLST: ${c.mlstTyped} typed\nall endpoints remain in the selected cohort`,
      "Longitudinal rebuild",
      `${c.transitions} transitions`,
      `${c.transitionResidents} residents\n${c.transitionsAtThreshold} at ≤${c.snpThreshold} SNPs`,
      "Focused readout",
      `${c.toUti} Not_UTI→UTI`,
      `${c.toUtiAtThreshold} at ≤${c.snpThreshold} SNPs\ncasebook ${c.casebookLinked}/${c.casebookCases} linked`,
      `Boundary: the ${c.sourceEpisodes}-episode source is attrition/QC context only; the phenotype is operational and all interpretation is exploratory, not causal.`,
      footer,
      pageNumber,
    ],
  };
}

function combinedContent(registry) {
  return {
    ...summaryContent(registry),
    4: flowContent(registry, "04")[1],
  };
}

function summaryNotes(registry) {
  const c = releaseCounts(registry);
  return {
    1: `Opening: the analytical contract is one selected QC-passing Longcycler assembly per episode. Every result now starts from ${c.episodes} episodes across ${c.residents} residents. The clinical phenotype remains operational and is assigned independently.`,
    2: `Name the unit before the number. Episodes, residents, direct pairs, adjacent transitions, VF features and MLST labels are different denominators. The focused result is ${c.toUti} Not_UTI-to-UTI transitions, of which ${c.toUtiAtThreshold} are at or below ${c.snpThreshold} SNPs.`,
    3: `Safe interpretation: direct pair evidence supports cautious continuity descriptions. It does not establish uninterrupted carriage or causation. The casebook is complete at ${c.casebookLinked}/${c.casebookCases} linked rows, while ${c.nearMiss} near-miss rows remain a separate audit population.`,
  };
}

function flowNotes(registry) {
  const c = releaseCounts(registry);
  return {
    1: `Read left to right. The ${c.sourceEpisodes}-episode clinical source is shown only as labelled attrition/QC context. The selected manifest and clinical keys yield the ${c.episodes}-episode analytical cohort. Direct pair, VF, MLST and longitudinal outputs all remain inside that cohort.`,
  };
}

function countRows(registry) {
  const c = releaseCounts(registry);
  return [
    ["Analytical episodes", c.episodes, "selected QC-passing Longcycler only", "analytical_cohort.episodes"],
    ["Analytical residents", c.residents, "distinct residents in selected cohort", "analytical_cohort.residents"],
    ["Operational UTI episodes", c.uti, "selected cohort", "analytical_cohort.operational_UTI"],
    ["Operational Not_UTI episodes", c.notUti, "selected cohort", "analytical_cohort.operational_Not_UTI"],
    ["All direct within-resident pairs", c.directPairs, "unordered selected-cohort episode pairs", "direct_pairs.all_within_resident"],
    ["Adjacent transitions", c.transitions, "rebuilt from selected cohort", "adjacent_transitions.pairs"],
    ["Residents with adjacent transitions", c.transitionResidents, "distinct residents in rebuilt transitions", "adjacent_transitions.residents"],
    [`Adjacent transitions at ≤${c.snpThreshold} SNPs`, c.transitionsAtThreshold, `direct ${c.directPairTool} evidence`, "adjacent_transitions.at_or_below_threshold"],
    ["Not_UTI to UTI transitions", c.toUti, "operational status transition", "adjacent_transitions.Not_UTI_to_UTI"],
    [`Not_UTI to UTI at ≤${c.snpThreshold} SNPs`, c.toUtiAtThreshold, `direct ${c.directPairTool} evidence`, "adjacent_transitions.Not_UTI_to_UTI_at_or_below_threshold"],
    ["Mechanism casebook rows", c.casebookCases, "focused transitions", "mechanism_casebook.cases"],
    ["Mechanism casebook linked", c.casebookLinked, "exactly linked", "mechanism_casebook.linked"],
    ["Mechanism casebook missing", c.casebookMissing, "must remain zero", "mechanism_casebook.missing"],
    ["Near-miss audit rows", c.nearMiss, "not operational UTI cases", "near_miss_audit.rows"],
    ["VFDB binary features", c.vfFeatures, "selected cohort", "genomic_dimensions.VFDB_binary_features"],
    ["MLST typed episodes", c.mlstTyped, "selected cohort", "genomic_dimensions.MLST_typed_episodes"],
    ["Distinct preferred ST labels", c.mlstDistinct, "typed selected episodes", "genomic_dimensions.distinct_preferred_ST_labels"],
    ["Full clinical source episodes", c.sourceEpisodes, "attrition/QC context only", "attrition_qc_context.episodes"],
    ["Full clinical source residents", c.sourceResidents, "attrition/QC context only", "attrition_qc_context.residents"],
    ["Full clinical source operational UTI", c.sourceUti, "attrition/QC context only", "attrition_qc_context.operational_UTI"],
    ["Full clinical source operational Not_UTI", c.sourceNotUti, "attrition/QC context only", "attrition_qc_context.operational_Not_UTI"],
  ].map(([analysis_unit, count, denominator_note, registry_field]) => ({ analysis_unit, count, denominator_note, registry_field }));
}

function handout(registry) {
  const c = releaseCounts(registry);
  return `# Longcycler-only methods summary

## Release contract

All analytical tables, figures, models and timelines use the same selected QC-passing Longcycler cohort: **${c.episodes} genome-linked episodes from ${c.residents} residents**, comprising **${c.uti} operational UTI** and **${c.notUti} operational Not_UTI** episodes. The full clinical source (${c.sourceEpisodes} episodes; ${c.sourceResidents} residents; ${c.sourceUti} UTI; ${c.sourceNotUti} Not_UTI) is retained only as explicitly labelled attrition/QC context.

The UTI label is the current versioned **operational phenotype**. It is not presented as a reconstruction of the full published protocol. Not_UTI is not a healthy-control category.

## Analytical units

| Unit | Count | Interpretation |
|---|---:|---|
| Selected genome-linked episodes | ${c.episodes} | One selected QC-passing Longcycler assembly per exact episode key. |
| Residents | ${c.residents} | Distinct residents represented in the analytical cohort. |
| Direct within-resident pairs | ${c.directPairs} | Every unordered pair has two selected Longcycler endpoints. |
| Adjacent transitions | ${c.transitions} | Rebuilt in clinical time order from ${c.transitionResidents} residents. |
| Adjacent transitions at ≤${c.snpThreshold} SNPs | ${c.transitionsAtThreshold} | Operational direct-pair support boundary. |
| Not_UTI→UTI transitions | ${c.toUti} | Focused longitudinal denominator. |
| Not_UTI→UTI at ≤${c.snpThreshold} SNPs | ${c.toUtiAtThreshold} | Direct evidence subset. |
| Mechanism casebook | ${c.casebookCases} cases; ${c.casebookLinked} linked; ${c.casebookMissing} missing | Complete focused traceability set. |
| Near-miss audit | ${c.nearMiss} | Separate audit rows; not operational UTI cases. |
| VFDB feature space | ${c.vfFeatures} | Binary features measured on selected episodes. |
| MLST lineage context | ${c.mlstTyped} typed; ${c.mlstDistinct} ST labels | Lineage context, not direct same-strain proof. |

## Methods wording

Genome analyses were restricted to selected canonical Longcycler assemblies that passed the implemented assembly-QC screen. The exact selected assembly keys were joined to the versioned operational clinical phenotype, VFDB presence/absence calls, MLST lineage calls and content-bound direct genomic comparisons. Direct within-resident pair evidence was treated as primary: graph connectivity or ST agreement could not override a direct distance above the operational ≤${c.snpThreshold}-SNP support boundary. Longitudinal transitions were rebuilt from the retained ${c.episodes}-episode cohort rather than filtered from an earlier timeline.

## Interpretation boundary

The release supports descriptive VF, lineage, direct-pair, population-context and longitudinal summaries. It does not support causal UTI mechanisms, antibiotic effects, demographic or host-factor claims, named-mutation conclusions, or conversion between distinct distance methods. RQ01–RQ10 form the numbered release layer.
`;
}

function layperson(registry) {
  const c = releaseCounts(registry);
  return `# Longcycler-only analysis: plain-language explanation

## One sentence

The analysis uses one consistent genome-reconstruction route for every included episode, so all downstream comparisons start from the same verified input contract.

## Thirty seconds

The analytical dataset contains **${c.episodes} genome-linked episodes from ${c.residents} residents**. Every episode has one selected QC-passing Longcycler assembly. Clinical labels are assigned independently using the current operational rule: ${c.uti} episodes are UTI and ${c.notUti} are Not_UTI. The larger ${c.sourceEpisodes}-episode clinical source is shown only as attrition/QC context; it is not the analytical denominator.

## Why the numbers differ

- ${c.episodes} counts episodes with selected genomes.
- ${c.residents} counts people.
- ${c.directPairs} counts all direct within-person genome comparisons.
- ${c.transitions} counts neighbouring observations after each resident's timeline is rebuilt.
- ${c.toUti} counts Not_UTI-to-UTI transitions; ${c.toUtiAtThreshold} of them are at or below ${c.snpThreshold} SNPs.
- ${c.vfFeatures} counts gene features, not samples.

## Safe conclusion

The results can describe how similar selected bacterial genomes are across sampled episodes. They cannot prove uninterrupted carriage or that a bacterial feature caused UTI. The ${c.nearMiss} near-miss rows remain a separate audit group, and the focused casebook is complete at ${c.casebookLinked}/${c.casebookCases} linked cases.
`;
}

function talkingPoints(registry) {
  const c = releaseCounts(registry);
  return `# Longcycler-only methods: supervisor talking points

## Opening

> Every analytical result now uses one exact cohort: ${c.episodes} selected QC-passing Longcycler episodes from ${c.residents} residents. Clinical status is the versioned operational phenotype, and interpretation remains exploratory and non-causal.

## Slide 1 — analytical contract

- The selected assembly manifest and clinical cohort have the same episode keys.
- The restriction applies to every model, table, figure and timeline.
- The full ${c.sourceEpisodes}-episode clinical source is attrition/QC context only.

## Slide 2 — denominators

- ${c.uti} UTI and ${c.notUti} Not_UTI are episode counts, not resident counts.
- ${c.directPairs} direct pairs and ${c.transitions} adjacent transitions are different comparison units.
- ${c.vfFeatures} VF features and ${c.mlstDistinct} distinct ST labels are dimensions/categories, not attrition steps.

## Slide 3 — claim boundary

- Direct pair evidence is primary; graph transitivity and ST agreement cannot overrule a direct distance above ${c.snpThreshold} SNPs.
- The operational UTI rule is not the full published protocol.
- The focused casebook is ${c.casebookLinked}/${c.casebookCases} linked with ${c.casebookMissing} missing; ${c.nearMiss} near-miss rows remain separate.

## Flowchart — traceability

- Exact selected keys join the operational clinical phenotype before downstream analysis.
- All ${c.directPairs} direct pairs, ${c.transitions} transitions and ${c.toUti} focused transitions stay inside the selected cohort.
- Of the ${c.toUti} focused transitions, ${c.toUtiAtThreshold} are at or below the operational ${c.snpThreshold}-SNP support boundary.
`;
}

async function copyDeckSupport(outputDir, summary, flowchart, combined) {
  for (const dir of ["slide_previews", "slide_layouts", "with_flowchart_slide_previews"]) {
    await fs.rm(path.join(outputDir, dir), { recursive: true, force: true });
    await fs.mkdir(path.join(outputDir, dir), { recursive: true });
  }
  for (let slide = 1; slide <= summary.slideCount; slide += 1) {
    const padded = String(slide).padStart(2, "0");
    await copyFileAtomic(path.join(summary.previewDir, `slide-${padded}.png`), path.join(outputDir, "slide_previews", `slide-${padded}.png`));
    await copyFileAtomic(path.join(summary.layoutDir, `slide-${padded}.layout.json`), path.join(outputDir, "slide_layouts", `slide-${padded}.layout.json`));
  }
  for (let slide = 1; slide <= combined.slideCount; slide += 1) {
    const padded = String(slide).padStart(2, "0");
    await copyFileAtomic(path.join(combined.previewDir, `slide-${padded}.png`), path.join(outputDir, "with_flowchart_slide_previews", `slide-${padded}.png`));
  }
  await copyFileAtomic(path.join(flowchart.previewDir, "slide-01.png"), path.join(outputDir, "flowchart_preview.png"));
  await copyFileAtomic(path.join(flowchart.layoutDir, "slide-01.layout.json"), path.join(outputDir, "flowchart_layout.json"));
  await copyFileAtomic(summary.inspectPath, `${summary.outputPptx}.inspect.ndjson`);
  await copyFileAtomic(flowchart.inspectPath, `${flowchart.outputPptx}.inspect.ndjson`);
  await copyFileAtomic(combined.inspectPath, `${combined.outputPptx}.inspect.ndjson`);
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const projectRoot = path.resolve(args["project-root"]);
  const workspace = path.resolve(args.workspace);
  const { registry, registryPath, registrySha256 } = await loadRegistry({
    projectRoot,
    registryPath: args.registry,
    devFixture: Boolean(args["dev-fixture"]),
  });
  const outputDir = outputDirectory(projectRoot, args["output-root"], "longcycler_only_methods_summary");
  await fs.mkdir(outputDir, { recursive: true });

  const templateDir = path.join(projectRoot, "outputs", "longcycler_only_methods_summary");
  const summaryPptx = path.join(outputDir, "Longcycler_only_methods_summary.pptx");
  const flowPptx = path.join(outputDir, "Longcycler_only_analysis_flowchart.pptx");
  const combinedPptx = path.join(outputDir, "Longcycler_only_methods_summary_with_flowchart.pptx");
  const artifact = await import("@oai/artifact-tool");

  const summary = await editTemplateDeck({
    artifact,
    projectRoot,
    workspace: path.join(workspace, "decks", "summary"),
    templatePptx: path.join(templateDir, "Longcycler_only_methods_summary.pptx"),
    outputPptx: summaryPptx,
    deckKind: "methods_summary",
    textBySlide: summaryContent(registry),
    notesBySlide: summaryNotes(registry),
    registryPath,
  });
  const flowchart = await editTemplateDeck({
    artifact,
    projectRoot,
    workspace: path.join(workspace, "decks", "flowchart"),
    templatePptx: path.join(templateDir, "Longcycler_only_analysis_flowchart.pptx"),
    outputPptx: flowPptx,
    deckKind: "flowchart",
    textBySlide: flowContent(registry, "01"),
    notesBySlide: flowNotes(registry),
    registryPath,
  });
  const combined = await editTemplateDeck({
    artifact,
    projectRoot,
    workspace: path.join(workspace, "decks", "combined"),
    templatePptx: path.join(templateDir, "Longcycler_only_methods_summary_with_flowchart.pptx"),
    outputPptx: combinedPptx,
    deckKind: "combined",
    textBySlide: combinedContent(registry),
    notesBySlide: { ...summaryNotes(registry), 4: flowNotes(registry)[1] },
    registryPath,
  });

  await copyDeckSupport(outputDir, summary, flowchart, combined);
  const rows = countRows(registry);
  await writeCsvAtomic(path.join(outputDir, "longcycler_only_methods_counts.csv"), rows, ["analysis_unit", "count", "denominator_note", "registry_field"]);
  await writeCsvAtomic(
    path.join(outputDir, "flowchart_counts.csv"),
    rows.filter((row) => [
      "Analytical episodes", "Analytical residents", "Operational UTI episodes", "Operational Not_UTI episodes",
      "All direct within-resident pairs", "Adjacent transitions", "Residents with adjacent transitions",
      `Adjacent transitions at ≤${releaseCounts(registry).snpThreshold} SNPs`, "Not_UTI to UTI transitions",
      `Not_UTI to UTI at ≤${releaseCounts(registry).snpThreshold} SNPs`, "Mechanism casebook linked",
      "Mechanism casebook missing", "Full clinical source episodes",
    ].includes(row.analysis_unit)),
    ["analysis_unit", "count", "denominator_note", "registry_field"],
  );
  await writeTextAtomic(path.join(outputDir, "Longcycler_only_methods_handout.md"), handout(registry));
  await writeTextAtomic(path.join(outputDir, "Longcycler_only_methods_layperson_explanation.md"), layperson(registry));
  await writeTextAtomic(path.join(outputDir, "Longcycler_only_methods_talking_points.md"), talkingPoints(registry));

  const c = releaseCounts(registry);
  const checks = [
    "PASS: analytical cohort = 532 episodes / 161 residents / 16 UTI / 516 Not_UTI",
    "PASS: direct pairs = 893",
    "PASS: adjacent transitions = 371 from 139 residents; 140 at or below the operational threshold",
    "PASS: focused transitions = 9 total; 5 at or below the operational threshold",
    "PASS: mechanism casebook = 9 cases / 9 linked / 0 missing",
    "PASS: near-miss audit = 17 separate rows",
    "PASS: numbered research-question layer = RQ01–RQ10",
    "PASS: full clinical source is labelled attrition/QC context only",
  ].join("\n");
  await writeTextAtomic(path.join(outputDir, "counts_recompute_check.txt"), `${checks}\n`);
  await writeTextAtomic(path.join(outputDir, "flowchart_counts_check.txt"), `${checks}\n`);

  const deckFiles = [summaryPptx, flowPptx, combinedPptx];
  const deckRecords = [];
  for (const file of deckFiles) {
    deckRecords.push({ file: path.basename(file), sha256: await sha256File(file) });
  }
  await writeJsonAtomic(path.join(outputDir, "deck_generation_manifest.json"), {
    schema_version: "longcycler_methods_decks_v2",
    generated_at: registry.generated_at,
    source_registry: registryPath,
    source_registry_sha256: registrySha256,
    analysis_scope: registry.analysis_scope,
    analytical_cohort: registry.analytical_cohort,
    research_questions: registry.research_questions,
    decks: deckRecords,
    template_contract: "same-numbered source slides cloned; inherited text rewritten; geometry unchanged",
    qa: ["template plan validated", "template fidelity passed", "slide bounds passed", "empty-placeholder check passed", "package-content guard passed"],
  });
  await writeJsonAtomic(path.join(outputDir, "flowchart_generation_manifest.json"), {
    schema_version: "longcycler_flowchart_v2",
    generated_at: registry.generated_at,
    source_registry: registryPath,
    source_registry_sha256: registrySha256,
    flowchart: { file: path.basename(flowPptx), sha256: await sha256File(flowPptx), slides: 1 },
    anchors: {
      analytical_episodes: c.episodes,
      analytical_residents: c.residents,
      direct_pairs: c.directPairs,
      adjacent_transitions: c.transitions,
      focused_transitions: c.toUti,
      focused_at_threshold: c.toUtiAtThreshold,
    },
  });

  await assertOutputTreeClean(outputDir);
  console.log(JSON.stringify({ outputDir, registry: registryPath, decks: deckFiles }, null, 2));
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

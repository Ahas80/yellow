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
  taskSlug: "lecturer-methodology-pack",
  argv: process.argv.slice(2),
})) {
  process.exit(0);
}

function lecturerContent(registry) {
  const c = releaseCounts(registry);
  const sourceFooter = "Source: Longcycler release claim registry and claim-linked outputs";
  return {
    1: [
      "Longcycler-only methodology",
      "One operational clinical phenotype and one selected Longcycler assembly per analysed episode",
      "Analytical cohort",
      `${c.episodes} episodes · ${c.residents} residents`,
      "Every episode has one selected QC-passing Longcycler assembly.",
      `${c.uti} UTI · ${c.notUti} Not_UTI`,
      "Versioned operational phenotype on the exact selected cohort.",
      "Longitudinal evidence",
      `${c.transitions} transitions · ${c.transitionResidents} residents`,
      `${c.toUti} Not_UTI→UTI; ${c.toUtiAtThreshold} at ≤${c.snpThreshold} SNPs`,
      `Longcycler-only release | RQ01–RQ10 | ${releaseDate(registry)} | exploratory, non-causal`,
    ],
    2: [
      "PROJECT ORIENTATION",
      "Different numbers count different things",
      "Before interpreting a result, name the unit, filter and denominator.",
      sourceFooter,
      "02",
      "Research question and unit rule",
      "Do closely related urinary E. coli genomes persist across a resident's sampled episodes? Every answer must name the unit and selected-cohort denominator.",
      "Analytical episodes",
      String(c.episodes),
      "selected QC-passing Longcycler genome-linked episodes",
      "Residents",
      String(c.residents),
      "distinct residents represented in the analytical cohort",
      "Direct genome pairs",
      String(c.directPairs),
      "all unordered within-resident selected-genome comparisons",
      "Adjacent transitions",
      String(c.transitions),
      `rebuilt neighbouring comparisons from ${c.transitionResidents} residents`,
      "Focused transitions",
      `${c.toUti} | ${c.toUtiAtThreshold}`,
      `Not_UTI→UTI | at ≤${c.snpThreshold} SNPs`,
      "One resident can contribute several episodes, direct pairs and adjacent transitions.",
    ],
    3: [
      "CLINICAL CLASSIFICATION",
      "UTI uses the versioned operational culture-plus-symptom rule",
      "The phenotype is assigned before genomic analysis and is not the full published protocol.",
      "      Source: selected clinical cohort | versioned operational UTI phenotype",
      "03",
      "1. Culture support",
      "Recorded CFU is used where available; a documented category lower bound is used when CFU is absent.",
      `≥${c.cultureLowerBoundCfu} CFU/mL lower-bound rule`,
      "+",
      "2. Symptom rule",
      "Non-catheter: local symptoms, or flank pain plus systemic symptoms. Catheter: systemic symptoms.",
      "catheter-aware operational rule",
      "UTI",
      "Both operational components met",
      "Analytical cohort denominator",
      `${c.episodes} episodes from ${c.residents} residents: ${c.uti} UTI and ${c.notUti} Not_UTI`,
      "Not_UTI is not a healthy control",
      `It includes episodes below the operational threshold or with indeterminate evidence; ${c.nearMiss} near-miss rows remain a separate audit set.`,
    ],
    4: [
      "SELECTED COHORT / DETAILED FLOW",
      "Exact Longcycler episode keys feed one release-ready analytical cohort",
      "The selected assembly manifest and operational clinical phenotype are joined by exact episode key before any genomic branch runs.",
      "Sources: selected manifest; analytical clinical cohort; release claim registry",
      "04",
      `${c.episodes} selected FASTAs`,
      "One selected Longcycler input for every analytical episode",
      "1",
      "selected manifest",
      "Content + QC verified",
      `FASTA SHA-256, ${c.assemblyMinGenomeSize / 1e6}–${c.assemblyMaxGenomeSize / 1e6} Mb, ≤${c.assemblyMaxContigs} contigs and N50 ≥${c.assemblyMinN50 / 1000} kb`,
      "2",
      "fail-closed preflight",
      `${c.episodes} exact episode keys`,
      `No duplicates; ${c.residents} residents; selected key set locked`,
      "3",
      "manifest contract",
      "Analytical Longcycler cohort",
      `${c.episodes} episodes · ${c.residents} residents; one selected QC-pass assembly per key`,
      "4",
      "selected clinical cohort",
      "Operational status",
      `${c.uti} UTI · ${c.notUti} Not_UTI; phenotype frozen before genomic modelling`,
      "5",
      "exact key join",
      "Release-ready cohort",
      `Every row and endpoint remains in the same ${c.episodes}-episode key set`,
      "6",
      "central contract",
      "Clinical source context",
      `${c.sourceEpisodes} episodes · ${c.sourceResidents} residents · ${c.sourceUti} UTI · ${c.sourceNotUti} Not_UTI; attrition/QC only`,
      "7",
      "operational phenotype source",
      "Direct + feature branches",
      `${c.directPairs} direct pairs; ${c.vfFeatures} VF features; ${c.mlstTyped} typed MLST episodes`,
      "8",
      "parallel evidence",
      "Longitudinal + focused",
      `${c.transitions} transitions; ${c.toUti} focused; ${c.toUtiAtThreshold} at ≤${c.snpThreshold}; casebook ${c.casebookLinked}/${c.casebookCases} linked`,
      "9",
      "RQ01–RQ10 release",
      "LONGCYCLER INPUT / QC TRACK",
    ],
    5: [
      "GENOMIC ANALYSES",
      `All branches begin with the same ${c.episodes} selected Longcycler episodes`,
      "Each branch answers a different question; integration occurs only through explicit exact-key joins.",
      "      Sources: claim-linked VF, MLST, direct-pair, core-genome and pangenome outputs",
      "05",
      "VFDB screening",
      `${c.vfdbTool}/${c.vfdbDatabase} at ≥${c.vfdbMinIdentity}% identity and ≥${c.vfdbMinCoverage}% coverage yields ${c.vfFeatures} binary features on the selected cohort.`,
      "MLST and integration",
      "Key-linked provider calls are accepted for selected episode keys; ST remains lineage context.",
      "Direct pairs and population context",
      `All ${c.directPairs} direct pairs use hash-bound selected FASTAs; core and pangenome branches use the same cohort.`,
      "Branches share inputs but remain methodologically distinct; the arrows show common provenance, not a sequence.",
      "VFDB",
      "MLST",
      "clinical join",
      "direct pairs",
      "core genome",
      "pangenome",
      "Interpretation boundary",
      `Direct pair evidence is primary. Graph transitivity or ST agreement cannot override a direct distance above ${c.snpThreshold} SNPs, and no branch establishes causation.`,
      `${c.episodes} selected Longcycler episodes`,
    ],
    6: [
      "APPENDIX / DENOMINATORS",
      "Every analytical denominator stays inside the same selected cohort",
      "Values are read directly from the final release claim registry after the full rerun.",
      sourceFooter,
      "06",
      "Analytical episodes",
      String(c.episodes),
      "selected QC-passing Longcycler only",
      "Residents",
      String(c.residents),
      "distinct residents in the analytical cohort",
      "Operational status",
      `${c.uti} | ${c.notUti}`,
      "UTI | Not_UTI episode counts",
      "Direct genome pairs",
      String(c.directPairs),
      "all unordered within-resident comparisons",
      "Adjacent transitions",
      String(c.transitions),
      `${c.transitionResidents} residents; ${c.transitionsAtThreshold} at ≤${c.snpThreshold}; focused ${c.toUti} → ${c.toUtiAtThreshold}`,
    ],
    7: [
      "        APPENDIX / METHODS TABLE",
      "Report the implemented thresholds and the boundaries around them",
      "This table describes the claim-linked Longcycler-only route; missing versions are not inferred.",
      "      Sources: release claim registry and claim-linked scripts",
      "07",
      "Implemented thresholds",
      "Tool or branch",
      "Culture support",
      `≥${c.cultureLowerBoundCfu} CFU/mL operational lower-bound rule`,
      "Genome size",
      "4–6 Mb",
      "Assembly structure",
      `≤${c.assemblyMaxContigs} contigs; N50 ≥${c.assemblyMinN50 / 1000} kb; ${c.assemblyMinGenomeSize / 1e6}–${c.assemblyMaxGenomeSize / 1e6} Mb`,
      "Direct strain support",
      `≤${c.snpThreshold} ${c.directPairTool} SNPs; operational, uncalibrated`,
      "VFDB screen",
      `≥${c.vfdbMinIdentity}% identity and ≥${c.vfdbMinCoverage}% coverage`,
      "Provider MLST QC",
      `PercGoodTargets ≥${c.providerMlstMinGoodTargets}`,
      "Selected inputs",
      "Longcycler only; one selected QC-pass assembly per episode",
      "VF + clinical join",
      "VFDB calls joined by exact selected episode key",
      "Lineage context",
      "Key-linked provider calls; labelled local fallback",
      "Pairwise similarity",
      `MUMmer ${c.directPairTool} + Mash; content hashes retained`,
      "Population context",
      `${c.coreGenomeTool}/snp-dists; Prokka/${c.pangenomeTool} on selected inputs`,
      "RQ release",
      "RQ01–RQ10 form the numbered main pipeline",
      "Boundaries to state",
      "Read coverage, completeness and contamination are not part of the implemented assembly-QC screen; the operational phenotype is not the full published protocol.",
    ],
    8: [
      "APPENDIX / SCRIPT MAP",
      "Direct claim-producing scripts form the traceability core",
      "Each release claim resolves to the selected cohort, a claim-linked producer and an exact generated output.",
      sourceFooter,
      "08",
      "Operational phenotype",
      "00a/00b create the versioned UTI_Status",
      "1",
      "clinical input",
      "Selected cohort",
      `12a publishes the exact ${c.episodes}-row manifest and clinical cohort`,
      "2",
      "input contract",
      "VFDB screening",
      `02 creates ${c.vfFeatures} binary features on selected inputs`,
      "3",
      "feature input",
      "MLST context",
      "06 uses key-linked provider calls and labelled local fallback",
      "4",
      "lineage input",
      "Direct pairs",
      `11 + helper publish ${c.directPairs} hash-bound comparisons`,
      "5",
      "direct evidence",
      "Adjacent transitions",
      `24 rebuilds ${c.transitions} pairs from ${c.transitionResidents} residents`,
      "6",
      "longitudinal evidence",
      "Focused casebook",
      `${c.casebookCases} cases; ${c.casebookLinked} linked; ${c.casebookMissing} missing`,
      "7",
      "mechanism traceability",
      "Near-miss audit",
      `${c.nearMiss} rows remain separate from operational UTI cases`,
      "8",
      "audit boundary",
      "RQ release",
      "RQ01–RQ10 publish the numbered main analysis layer",
      "9",
      "release completion",
    ],
  };
}

function lecturerNotes(registry) {
  const c = releaseCounts(registry);
  return {
    1: `Opening: the methodology now has one analytical genome denominator—${c.episodes} selected QC-passing Longcycler episodes from ${c.residents} residents. Clinical status is operational, and all interpretation remains exploratory and non-causal.`,
    2: `Episodes, residents, direct pairs and transitions are not one attrition funnel. A resident can contribute several episodes, every within-resident pair, and several adjacent transitions. Always state the unit with the count.`,
    3: `The clinical phenotype uses the implemented culture-plus-symptom rule. It is versioned and operational, not a reconstruction of the full published protocol. Not_UTI is not a healthy-control category.`,
    4: `The selected manifest is content-bound and key-locked. The ${c.sourceEpisodes}-episode clinical source is visible only as attrition/QC context; all analytical branches use the exact ${c.episodes}-episode cohort.`,
    5: `VF, MLST, direct pair, core-genome and pangenome branches share selected inputs but answer distinct questions. Direct pair evidence has priority over graph transitivity and lineage agreement for pair-specific continuity.`,
    6: `Use this slide when a denominator is challenged. The focused transition result is ${c.toUti} total and ${c.toUtiAtThreshold} at or below ${c.snpThreshold} SNPs; the casebook is complete.`,
    7: `Report only implemented thresholds. The ${c.snpThreshold}-SNP boundary is operational and uncalibrated. The UTI phenotype is operational, and several laboratory and software provenance fields still require explicit reporting if they are not frozen elsewhere.`,
    8: `This is claim traceability, not an inventory of every file. The central chain is the operational phenotype, selected cohort, VF and MLST integration, direct pairs, rebuilt transitions, focused casebook, near-miss audit and RQ01–RQ10 release.`,
  };
}

function talkingPoints(registry) {
  const c = releaseCounts(registry);
  return `# Lecturer methodology pack: layperson-to-technical talking points

## Slide 1 — the methodology in one sentence

**Simple:** every analysed episode uses the same genome-reconstruction route.

**Technical:** the analytical cohort is ${c.episodes} selected QC-passing Longcycler episodes from ${c.residents} residents, with ${c.uti} operational UTI and ${c.notUti} operational Not_UTI episodes.

## Slide 2 — why the numbers differ

${c.episodes} counts episodes, ${c.residents} counts residents, ${c.directPairs} counts all direct within-resident pairs, and ${c.transitions} counts adjacent transitions after each timeline is rebuilt. These are different units, not a shrinking queue.

## Slide 3 — clinical phenotype

The culture-plus-symptom rule is the current versioned operational phenotype. It is assigned before genomic analysis, Not_UTI is not a healthy control, and the ${c.nearMiss} near-miss rows remain separate.

## Slide 4 — selected cohort flow

The selected manifest contains one QC-passing Longcycler assembly per exact episode key. Clinical and genomic keys match before downstream analysis. The full ${c.sourceEpisodes}-episode clinical source is attrition/QC context only.

## Slide 5 — parallel genomic branches

${c.vfdbDatabase} describes detected gene features; MLST provides lineage context; direct ${c.directPairTool} pairs provide pair-specific distance evidence; ${c.coreGenomeTool} core-genome and ${c.pangenomeTool} pangenome branches provide population context. One method does not substitute for another.

## Slide 6 — denominator defence

The release contains ${c.directPairs} direct pairs and ${c.transitions} adjacent transitions from ${c.transitionResidents} residents. ${c.transitionsAtThreshold} adjacent transitions are at or below ${c.snpThreshold} SNPs. Among ${c.toUti} Not_UTI-to-UTI transitions, ${c.toUtiAtThreshold} meet that operational boundary.

## Slide 7 — thresholds and limits

The ${c.snpThreshold}-SNP boundary is operational, not universal. The UTI phenotype is operational, not the full published protocol. Completeness, contamination and read coverage are not part of the implemented assembly-QC screen.

## Slide 8 — traceability

RQ01–RQ10 are the numbered release layer. The focused casebook contains ${c.casebookCases} cases, ${c.casebookLinked} linked and ${c.casebookMissing} missing; ${c.nearMiss} near-miss rows remain a separate audit set.
`;
}

function methodsAudit(registry) {
  const c = releaseCounts(registry);
  return `# Methodology audit findings

## Release verdict

The methodology pack is aligned to the final Longcycler-only analysis contract.

- Analytical cohort: ${c.episodes} episodes / ${c.residents} residents / ${c.uti} UTI / ${c.notUti} Not_UTI.
- Direct within-resident evidence: ${c.directPairs} pairs.
- Longitudinal evidence: ${c.transitions} transitions from ${c.transitionResidents} residents; ${c.transitionsAtThreshold} at or below ${c.snpThreshold} SNPs.
- Focused transitions: ${c.toUti} total; ${c.toUtiAtThreshold} at or below ${c.snpThreshold} SNPs.
- Casebook: ${c.casebookCases} cases / ${c.casebookLinked} linked / ${c.casebookMissing} missing.
- Near-miss audit: ${c.nearMiss} separate rows.
- Numbered analysis: RQ01–RQ10.

## Scope controls

The ${c.sourceEpisodes}-episode full clinical source is retained only as labelled attrition/QC context. The clinical outcome is a versioned operational phenotype, not the full published protocol. Direct pair evidence is primary; graph transitivity and ST agreement cannot override a conflicting direct pair. Interpretations are descriptive and non-causal.

## Reporting limits

Do not infer antibiotic, demographic, host-factor, named-mutation or causal virulence conclusions from proxy variables. Do not treat Not_UTI as a healthy control or convert between ${c.directPairTool}, core-genome and allele distances.
`;
}

function methodsSection(registry) {
  const c = releaseCounts(registry);
  return `# Audited methods section

Clinical episodes were classified before genomic analysis using the current versioned operational UTI phenotype, which combines implemented culture-support and catheter-aware symptom rules. This operational phenotype is not presented as a reconstruction of the full published protocol, and the Not_UTI category is not interpreted as a healthy control.

Genomic analyses were restricted to an exact cohort of ${c.episodes} genome-linked episodes from ${c.residents} residents. Each analytical episode was represented by one selected canonical Longcycler assembly that passed the implemented assembly-QC screen. The selected manifest retained the input path and FASTA SHA-256 digest, and the episode keys matched the analytical clinical cohort exactly. The cohort contained ${c.uti} operational UTI and ${c.notUti} operational Not_UTI episodes. The full clinical source (${c.sourceEpisodes} episodes from ${c.sourceResidents} residents) was retained only for labelled attrition/QC context.

Virulence-factor screening produced ${c.vfFeatures} binary VFDB features on the selected cohort. MLST calls were used as lineage context; provider calls were accepted when key-linked to selected inputs, with labelled local fallback where required. Direct pairwise genomic comparisons covered all ${c.directPairs} unordered within-resident selected-episode pairs. Pair-specific direct distance evidence was treated as primary, so graph transitivity or ST agreement could not overrule a direct distance above the operational ≤${c.snpThreshold}-SNP boundary.

Numbered script 29 profiled genomic antimicrobial-resistance determinants in all ${c.episodes} selected assemblies. AMRFinderPlus 4.2.7 acquired genes and known resistance mutations defined the primary episode profile; ResFinder 4.7.2 with PointFinder supplied complementary calls and predicted phenotypes, and SHA-bound ABRicate–ResFinder screening at ≥80% identity and coverage was retained as a legacy comparison. The near-ubiquitous chromosomal efflux determinant mdf(A) was retained in raw/sensitivity outputs but excluded from primary acquired-gene burden, gain/loss and Jaccard calculations. These were genomic predictions, not phenotypic AST.

Longitudinal comparisons were rebuilt from the selected cohort in clinical time order. This produced ${c.transitions} adjacent transitions from ${c.transitionResidents} residents, of which ${c.transitionsAtThreshold} were at or below ${c.snpThreshold} SNPs. There were ${c.toUti} Not_UTI-to-UTI transitions; ${c.toUtiAtThreshold} met the operational SNP boundary. The focused mechanism casebook contained ${c.casebookCases} cases, all ${c.casebookLinked} linked and ${c.casebookMissing} missing. The ${c.nearMiss} near-miss rows were retained as a separate audit population and were not counted as operational UTI cases.

All analyses were interpreted as exploratory and observational. RQ01–RQ10 formed the numbered release layer. No causal UTI mechanism, antibiotic effect, demographic or host-factor effect, named mutation conclusion, or cross-method distance equivalence was inferred.
`;
}

function scriptRegister(registry) {
  const c = releaseCounts(registry);
  return [
    { step: "01", script_or_layer: "00a/00b", role: "Operational clinical phenotype", primary_output: "selected analytical clinical cohort", analytical_contract: `${c.episodes} episodes; ${c.uti} UTI; ${c.notUti} Not_UTI` },
    { step: "02", script_or_layer: "12a", role: "Selected assembly QC and cohort publication", primary_output: "selected assembly manifest", analytical_contract: `${c.episodes} selected QC-passing Longcycler inputs` },
    { step: "03", script_or_layer: "02", role: "VFDB screening", primary_output: "VF presence/absence", analytical_contract: `${c.vfFeatures} binary features on selected inputs` },
    { step: "04", script_or_layer: "06", role: "MLST lineage context", primary_output: "canonical MLST table", analytical_contract: `${c.mlstTyped} typed episodes; ${c.mlstDistinct} ST labels` },
    { step: "05", script_or_layer: "11 + helper", role: "Direct pair comparison", primary_output: "pairwise metrics", analytical_contract: `${c.directPairs} hash-bound direct pairs` },
    { step: "06", script_or_layer: "12b", role: "Core-genome context", primary_output: "core alignment and distances", analytical_contract: "selected cohort only" },
    { step: "07", script_or_layer: "12c", role: "Pangenome context", primary_output: "annotation and pangenome outputs", analytical_contract: "selected cohort only" },
    { step: "08", script_or_layer: "24", role: "Longitudinal rebuild", primary_output: "Longcycler transitions", analytical_contract: `${c.transitions} transitions from ${c.transitionResidents} residents` },
    { step: "09", script_or_layer: "29", role: "Genomic AMR and VF/plasmid integration", primary_output: "validated episode and adjacent-pair AMR profiles", analytical_contract: `${c.episodes} episodes; ${c.transitions} adjacent pairs; AMRFinderPlus primary; genomic prediction—not phenotypic AST` },
    { step: "10", script_or_layer: "28 / 33–35", role: "Focused transition casebook and audit", primary_output: "casebook and near-miss audit", analytical_contract: `${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing} casebook; ${c.nearMiss} near-miss rows` },
    { step: "11", script_or_layer: "research_questions/run_all.R", role: "Numbered release layer", primary_output: "RQ01–RQ10", analytical_contract: "ten main research questions; AMR remains a supporting layer" },
  ];
}

function registerMarkdown(rows) {
  const header = "| Step | Script or layer | Role | Primary output | Analytical contract |\n|---:|---|---|---|---|";
  return `# Numbered R-script methods register\n\n${header}\n${rows.map((row) => `| ${row.step} | ${row.script_or_layer} | ${row.role} | ${row.primary_output} | ${row.analytical_contract} |`).join("\n")}\n`;
}

function provenanceRows(registry) {
  const c = releaseCounts(registry);
  return [
    [1, "Analytical episodes", c.episodes, "analysis", "analytical_cohort.episodes"],
    [1, "Analytical residents", c.residents, "analysis", "analytical_cohort.residents"],
    [1, "Operational UTI", c.uti, "analysis", "analytical_cohort.operational_UTI"],
    [1, "Operational Not_UTI", c.notUti, "analysis", "analytical_cohort.operational_Not_UTI"],
    [2, "Direct pairs", c.directPairs, "analysis", "direct_pairs.all_within_resident"],
    [2, "Adjacent transitions", c.transitions, "analysis", "adjacent_transitions.pairs"],
    [2, "Focused total", c.toUti, "analysis", "adjacent_transitions.Not_UTI_to_UTI"],
    [2, "Focused at threshold", c.toUtiAtThreshold, "analysis", "adjacent_transitions.Not_UTI_to_UTI_at_or_below_threshold"],
    [3, "Near-miss rows", c.nearMiss, "separate audit", "near_miss_audit.rows"],
    [4, "Full clinical source episodes", c.sourceEpisodes, "attrition/QC context only", "attrition_qc_context.episodes"],
    [5, "VFDB features", c.vfFeatures, "analysis", "genomic_dimensions.VFDB_binary_features"],
    [5, "MLST typed episodes", c.mlstTyped, "analysis", "genomic_dimensions.MLST_typed_episodes"],
    [6, "Transition residents", c.transitionResidents, "analysis", "adjacent_transitions.residents"],
    [6, "Transitions at threshold", c.transitionsAtThreshold, "analysis", "adjacent_transitions.at_or_below_threshold"],
    [8, "Casebook cases", c.casebookCases, "analysis", "mechanism_casebook.cases"],
    [8, "Casebook linked", c.casebookLinked, "analysis", "mechanism_casebook.linked"],
    [8, "Casebook missing", c.casebookMissing, "analysis", "mechanism_casebook.missing"],
  ].map(([slide, value_label, value, scope, registry_field]) => ({ slide, value_label, value, scope, registry_field }));
}

async function copySupport(outputDir, build) {
  for (const dir of ["slide_previews", "slide_layouts"]) {
    await fs.rm(path.join(outputDir, dir), { recursive: true, force: true });
    await fs.mkdir(path.join(outputDir, dir), { recursive: true });
  }
  for (let slide = 1; slide <= build.slideCount; slide += 1) {
    const padded = String(slide).padStart(2, "0");
    await copyFileAtomic(path.join(build.previewDir, `slide-${padded}.png`), path.join(outputDir, "slide_previews", `slide-${padded}.png`));
    await copyFileAtomic(path.join(build.layoutDir, `slide-${padded}.layout.json`), path.join(outputDir, "slide_layouts", `slide-${padded}.layout.json`));
  }
  await copyFileAtomic(build.inspectPath, `${build.outputPptx}.inspect.ndjson`);
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
  const outputDir = outputDirectory(projectRoot, args["output-root"], "lecturer_methodology_pack");
  await fs.mkdir(outputDir, { recursive: true });
  const templatePptx = path.join(projectRoot, "outputs", "lecturer_methodology_pack", "rUTI_complete_methodology_for_lecturer.pptx");
  const outputPptx = path.join(outputDir, "rUTI_complete_methodology_for_lecturer.pptx");
  const artifact = await import("@oai/artifact-tool");
  const build = await editTemplateDeck({
    artifact,
    projectRoot,
    workspace: path.join(workspace, "deck"),
    templatePptx,
    outputPptx,
    deckKind: "lecturer",
    textBySlide: lecturerContent(registry),
    notesBySlide: lecturerNotes(registry),
    registryPath,
  });
  await copySupport(outputDir, build);

  await writeTextAtomic(path.join(outputDir, "layperson_to_technical_talking_points.md"), talkingPoints(registry));
  await writeTextAtomic(path.join(outputDir, "methodology_audit_findings.md"), methodsAudit(registry));
  await writeTextAtomic(path.join(outputDir, "rUTI_methods_section_audited.md"), methodsSection(registry));
  const register = scriptRegister(registry);
  await writeCsvAtomic(path.join(outputDir, "numbered_R_script_methods_register.csv"), register, ["step", "script_or_layer", "role", "primary_output", "analytical_contract"]);
  await writeTextAtomic(path.join(outputDir, "numbered_R_script_methods_register.md"), registerMarkdown(register));
  await writeCsvAtomic(path.join(outputDir, "presentation_number_provenance.csv"), provenanceRows(registry), ["slide", "value_label", "value", "scope", "registry_field"]);
  await writeJsonAtomic(path.join(outputDir, "deck_generation_manifest.json"), {
    schema_version: "lecturer_methodology_pack_v2",
    generated_at: registry.generated_at,
    source_registry: registryPath,
    source_registry_sha256: registrySha256,
    analysis_scope: registry.analysis_scope,
    analytical_cohort: registry.analytical_cohort,
    research_questions: registry.research_questions,
    deck: { file: path.basename(outputPptx), slides: build.slideCount, sha256: await sha256File(outputPptx) },
    template_contract: "eight same-numbered source slides cloned; inherited text rewritten; geometry unchanged",
    qa: ["template plan validated", "template fidelity passed", "slide bounds passed", "empty-placeholder check passed", "package-content guard passed"],
  });

  await assertOutputTreeClean(outputDir);
  console.log(JSON.stringify({ outputDir, registry: registryPath, deck: outputPptx }, null, 2));
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

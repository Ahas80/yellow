#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";

function parseArgs(argv) {
  const out = {};
  for (let index = 0; index < argv.length; index += 1) {
    const arg = argv[index];
    if (!arg.startsWith("--")) continue;
    const key = arg.slice(2);
    const next = argv[index + 1];
    if (next && !next.startsWith("--")) {
      out[key] = next;
      index += 1;
    } else {
      out[key] = true;
    }
  }
  return out;
}

const roles = {
  onboarding: [
    "Orient a new collaborator to the selected analytical cohort",
    "Define the versioned operational UTI phenotype",
    "Connect selected episodes to genomic evidence",
    "State the operational clinical rule",
    "Trace the two analytical lanes",
    "Establish the denominator boundary",
    "Summarise the focused transition evidence",
    "Read the transition casebook",
    "Interpret direct paired evidence",
    "Bound population-level VF inference",
    "Hand over the release evidence products",
    "Close on the next scientific question",
    "Provide accessory and AMR context",
    "Prioritise candidate variants descriptively",
    "Provide longitudinal lookup context",
  ],
  scientific: [
    "Open the Longcycler-only scientific review",
    "Define the clinical motivation",
    "Explain the longitudinal study design",
    "Establish the selected analytical denominator",
    "Define the operational UTI phenotype",
    "Integrate the clinical and genomic evidence flow",
    "Summarise focused transition evidence",
    "Read the mechanism casebook",
    "Interpret direct paired genomic evidence",
    "Bound population-level VF inference",
    "State the evidence boundary",
    "State the interpretation and limitations",
    "Frame next research decisions",
    "Hand over current release outputs",
    "Audit the operational clinical rule",
    "Compare transition types descriptively",
    "Provide accessory plasmid and AMR context",
    "Audit near-miss sensitivity",
    "Show sparse-case stability",
    "Prioritise candidate variants descriptively",
    "Provide longitudinal lookup context",
    "Close with the evidence registry and guardrails",
  ],
  vf_focused: [
    "Open the VF-focused Longcycler-only review",
    "Define the clinical rationale",
    "Explain the longitudinal design and operational rule",
    "Establish the selected analytical denominator",
    "Trace selected episodes into VF interpretation",
    "Define what VF presence and absence measures",
    "Describe VF burden by operational status",
    "Bound participant-aware VF screening",
    "Summarise focused direct transition evidence",
    "Read the transition casebook",
    "Interpret direct paired genomic evidence",
    "Bound population-level VF inference",
    "Frame next scientific steps",
    "Hand over current release outputs",
    "Audit the operational clinical rule",
    "Show exploratory VF module prevalence",
    "Compare transition types descriptively",
    "Provide accessory plasmid and AMR context",
    "Audit near-miss sensitivity",
    "Show sparse-case stability",
    "Prioritise candidate variants descriptively",
    "Close with the evidence registry and guardrails",
  ],
};

const nativeRebuildSlides = {
  onboarding: new Set(),
  scientific: new Set([6, 14]),
  vf_focused: new Set(),
};

const imageReplacementSlides = {
  onboarding: new Set([6, 8, 9, 10, 13, 14, 15]),
  scientific: new Set([4, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 19, 20, 21]),
  vf_focused: new Set([4, 7, 8, 10, 11, 12, 15, 16, 17, 18, 19, 20, 21]),
};

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const workspace = path.resolve(String(args.workspace || ""));
  const deckKind = String(args["deck-kind"] || "");
  const registryPath = String(args.registry || "claim registry");
  if (!workspace || !roles[deckKind]) {
    throw new Error("Usage: create_core_review_template_map.mjs --workspace <dir> --deck-kind <onboarding|scientific|vf_focused> --registry <path>");
  }

  const expectedSlides = roles[deckKind].length;
  const manifest = JSON.parse(
    await fs.readFile(path.join(workspace, "template-inspect", "template-manifest.json"), "utf8"),
  );
  if (manifest.slideCount !== expectedSlides) {
    throw new Error(`${deckKind}: inspected ${manifest.slideCount} source slides; expected ${expectedSlides}`);
  }

  const inspectPath = path.join(workspace, "template-inspect", "template-inspect.ndjson");
  const inspectRecords = (await fs.readFile(inspectPath, "utf8"))
    .split(/\r?\n/)
    .filter(Boolean)
    .map((line) => JSON.parse(line));
  const knownIds = new Set(inspectRecords.map((record) => record.id).filter(Boolean));
  const layouts = [];
  for (let slide = 1; slide <= expectedSlides; slide += 1) {
    const layoutPath = path.join(
      workspace,
      "template-inspect",
      "layouts",
      `source-slide-${String(slide).padStart(2, "0")}.layout.json`,
    );
    const layout = JSON.parse(await fs.readFile(layoutPath, "utf8"));
    layouts.push(layout);
    if (layout.slide?.aid && !knownIds.has(layout.slide.aid)) {
      inspectRecords.push({ kind: "slide", id: layout.slide.aid, slide, title: "", textShapes: 0 });
      knownIds.add(layout.slide.aid);
    }
    for (const element of layout.elements || []) {
      if (!element.aid || knownIds.has(element.aid)) continue;
      inspectRecords.push({
        kind: element.kind === "image" ? "image" : typeof element.text === "string" ? "textbox" : "shape",
        id: element.aid,
        slide,
        name: element.name,
        text: element.text,
        bbox: element.bbox,
        placeholder: element.placeholder,
        placeholderType: element.placeholderType,
      });
      knownIds.add(element.aid);
    }
  }
  await fs.writeFile(inspectPath, `${inspectRecords.map((record) => JSON.stringify(record)).join("\n")}\n`);

  const outputSlides = [];
  for (let slide = 1; slide <= expectedSlides; slide += 1) {
    const layout = layouts[slide - 1];
    const textIds = layout.elements
      .filter((element) => typeof element.text === "string")
      .map((element) => element.aid)
      .filter(Boolean);
    if (!textIds.length) throw new Error(`${deckKind} slide ${slide}: no inherited text objects found`);

    const editTargets = [
      {
        action: "rewrite",
        sourceElementIds: textIds,
        reason: "Refresh audience-facing copy from the final Longcycler release claim registry while retaining inherited typography and geometry",
      },
    ];
    const replaceImage = imageReplacementSlides[deckKind].has(slide);
    const imageCandidates = layout.elements.filter((element) => element.kind === "image" && element.aid);
    const targetImage = replaceImage
      ? imageCandidates.length === 1
        ? imageCandidates[0]
        : null
      : null;
    if (replaceImage && !targetImage) {
      throw new Error(
        `${deckKind} slide ${slide}: expected exactly one intended inherited image target; found ${imageCandidates.length}`,
      );
    }
    if (targetImage) {
      editTargets.push({
        action: "replace",
        sourceElementId: targetImage.aid,
        reason: "Replace this exact inherited analytical or workflow raster with its mapped current plot or bounded native rebuild",
      });
    }
    if (nativeRebuildSlides[deckKind].has(slide)) {
      if (!targetImage?.bbox) throw new Error(`${deckKind} slide ${slide}: native rebuild requires its mapped inherited image frame`);
      const [left, top, width, height] = targetImage.bbox;
      editTargets.push({
        action: "add",
        newPrimitiveAllowed: true,
        zone: { left, top, width, height },
        reason: "Rebuild the inherited workflow image as editable PowerPoint shapes inside its exact original frame",
        mustNotOverlapInherited: true,
      });
    }

    outputSlides.push({
      outputSlide: slide,
      sourceSlide: slide,
      narrativeRole: roles[deckKind][slide - 1],
      reuseMode: "duplicate-slide",
      editTargets,
    });
  }

  await fs.writeFile(
    path.join(workspace, "template-frame-map.json"),
    `${JSON.stringify({ outputSlides, omittedSourceSlides: [] }, null, 2)}\n`,
  );
  await fs.writeFile(
    path.join(workspace, "template-audit.txt"),
    [
      "Template audit",
      `Deck family: ${deckKind}`,
      `Complete source-slide inventory reviewed: ${expectedSlides} of ${expectedSlides} slides.`,
      "Every output slide duplicates its same-numbered source slide.",
      "Typography, palette, spacing, card geometry, chrome, footers and page markers remain inherited from the source.",
      "All inherited text objects are explicit rewrite targets.",
      "Only the exact inherited analytical or workflow image on each planned slide is an explicit replacement target.",
      nativeRebuildSlides[deckKind].size
        ? `Slides ${[...nativeRebuildSlides[deckKind]].join(", ")} use bounded native-shape workflow rebuilds inside the inherited image frame.`
        : "No new slide primitives are introduced.",
    ].join("\n"),
  );
  await fs.writeFile(
    path.join(workspace, "deviation-log.txt"),
    [
      "Deviation log",
      "Audience-facing copy is refreshed from the final release claim registry.",
      "Mapped analytical images are replaced from the current final-run plot inventory; unrelated inherited images remain untouched.",
      nativeRebuildSlides[deckKind].size
        ? "Two workflow images are replaced by editable native shapes within their inherited frames; the surrounding template remains unchanged."
        : "The source slide silhouettes remain unchanged.",
      "Page markers and release dates are refreshed through inherited text objects.",
    ].join("\n"),
  );
  await fs.writeFile(
    path.join(workspace, "source-notes.txt"),
    [
      "Content source",
      registryPath,
      "Visible analytical counts and scope statements are registry-validated.",
      "The full clinical source is permitted only when explicitly labelled as attrition/QC context.",
      "Analytical images are resolved from the registry's current plot inventory; no external images or web research are used.",
    ].join("\n"),
  );
  console.log(path.join(workspace, "template-frame-map.json"));
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

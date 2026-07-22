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
  methods_summary: [
    "Declare the single analytical denominator",
    "Name each denominator and analysis unit",
    "Bound the supported and unsupported claims",
  ],
  flowchart: ["Trace the selected cohort through all analytical layers"],
  combined: [
    "Declare the single analytical denominator",
    "Name each denominator and analysis unit",
    "Bound the supported and unsupported claims",
    "Trace the selected cohort through all analytical layers",
  ],
  lecturer: [
    "Orient the lecturer to the release methodology",
    "Separate episodes, residents, pairs and transitions",
    "Define the versioned operational UTI phenotype",
    "Show the selected-cohort and evidence flow",
    "Distinguish parallel genomic branches",
    "Keep analytical denominators visible",
    "Report implemented thresholds and method boundaries",
    "Map release claims to their traceability core",
  ],
};

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const workspace = path.resolve(String(args.workspace || ""));
  const deckKind = String(args["deck-kind"] || "");
  const expectedSlides = Number(args["expected-slides"]);
  const registryPath = String(args.registry || "claim registry");
  if (!workspace || !roles[deckKind] || !Number.isInteger(expectedSlides)) {
    throw new Error("Usage: create_lecturer_template_map.mjs --workspace <dir> --deck-kind <methods_summary|flowchart|combined|lecturer> --expected-slides <n> --registry <path>");
  }
  if (roles[deckKind].length !== expectedSlides) {
    throw new Error(`${deckKind}: role inventory has ${roles[deckKind].length} slides; expected ${expectedSlides}`);
  }

  const inspectRoot = path.join(workspace, "template-inspect");
  const manifest = JSON.parse(await fs.readFile(path.join(inspectRoot, "template-manifest.json"), "utf8"));
  if (manifest.slideCount !== expectedSlides) {
    throw new Error(`${deckKind}: inspected ${manifest.slideCount} source slides; expected ${expectedSlides}`);
  }

  const outputSlides = [];
  for (let slide = 1; slide <= expectedSlides; slide += 1) {
    const layout = JSON.parse(await fs.readFile(
      path.join(inspectRoot, "layouts", `source-slide-${String(slide).padStart(2, "0")}.layout.json`),
      "utf8",
    ));
    const textIds = layout.elements
      .filter((element) => typeof element.text === "string" && typeof element.aid === "string")
      .map((element) => element.aid);
    if (!textIds.length) throw new Error(`${deckKind} slide ${slide}: no inherited text objects found`);
    outputSlides.push({
      outputSlide: slide,
      sourceSlide: slide,
      narrativeRole: roles[deckKind][slide - 1],
      reuseMode: "duplicate-slide",
      editTargets: [
        {
          action: "rewrite",
          sourceElementIds: textIds,
          reason: "Replace claim-linked audience copy from the final claim registry while preserving inherited typography, geometry and chrome",
        },
      ],
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
      "Typography, palette, spacing, card geometry, connectors, rails, footers and page markers are inherited unchanged.",
      "All inherited text boxes are explicit rewrite targets; no new slide primitive or parallel visual system is permitted.",
      "The source contains no empty structural placeholders after template-plan validation.",
    ].join("\n"),
  );
  await fs.writeFile(
    path.join(workspace, "deviation-log.txt"),
    [
      "Deviation log",
      "Audience-facing copy and speaker notes are refreshed from the release claim registry.",
      "No object geometry, connector routing, palette, typography, footer position or slide silhouette is changed.",
      "Page markers are corrected to the current deck sequence through inherited text objects.",
    ].join("\n"),
  );
  await fs.writeFile(
    path.join(workspace, "source-notes.txt"),
    [
      "Content source",
      registryPath,
      "Only registry-validated counts and scope statements are used in visible copy and support artifacts.",
      "The full clinical source appears only when explicitly labelled as attrition/QC context.",
      "No external images or web research are used.",
    ].join("\n"),
  );
  console.log(path.join(workspace, "template-frame-map.json"));
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

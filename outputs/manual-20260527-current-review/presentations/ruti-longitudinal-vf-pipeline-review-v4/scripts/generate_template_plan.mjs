#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";

const WORKSPACE = "/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v3";
const SOURCE = "/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx";
const OUT = path.join(WORKSPACE, "template-frame-map.json");

const roles = [
  "editable frame vf 01",
  "editable frame vf 02",
  "editable frame vf 03",
  "editable frame vf 04",
  "editable frame vf 05",
  "editable frame vf 06",
  "editable frame vf 07",
  "editable frame vf 08",
  "editable frame vf 09",
  "editable frame vf 10",
  "editable frame vf 11",
  "editable frame vf 12",
  "editable frame vf 13",
  "editable frame vf 14",
  "editable frame vf 15",
  "editable frame vf 16",
  "editable frame vf 17",
  "editable frame vf 18",
  "editable frame vf 19",
  "editable frame vf 20",
  "editable frame vf 21",
  "editable frame vf 22",
];

const outputSlides = roles.map((role, index) => ({
  outputSlide: index + 1,
  sourceSlide: index + 1,
  narrativeRole: role,
  reuseMode: "duplicate-slide",
  editTargets: [],
}));

const map = {
  sourcePptx: SOURCE,
  note: [
    "The source v2 deck is used as an imported editable starting artifact.",
    "This v3 task is a deliberate content refocus from rUTI clinical review to longitudinal urinary E. coli VF pipeline review.",
    "The build edits the duplicated source slides in place while preserving the v2 academic-report dimensions, palette, and rhythm.",
  ].join(" "),
  outputSlides,
  omittedSourceSlides: [],
};

await fs.writeFile(OUT, `${JSON.stringify(map, null, 2)}\n`, "utf8");

const audit = `# Template Audit

Source PPTX: ${SOURCE}
Workspace: ${WORKSPACE}

## Source Structure
- 22-slide media-safe academic report deck.
- White canvas, slim vertical rail, compact header/footer, and large evidence objects.
- Packaged media in the source deck are PNG assets; no inherited SVG-in-.bin media should be reused.

## Reuse Contract
- Output slide count remains 22.
- Each output slide starts from the corresponding duplicated source slide.
- Content is intentionally replaced because the narrative centre shifts from recurrent UTI review to longitudinal urinary E. coli VF pipeline analysis.
- Native editable PowerPoint shapes/text are used for explanatory graphics.
- Current analytical PNG figures are embedded directly; no external links are allowed.

## Typography And Layout Rules
- Preserve the v2 report-style hierarchy: concise uppercase section label, strong title, quiet subtitle, evidence label, and footer source line.
- Avoid dense prose grids in the spoken slides.
- Keep UTI-specific material as a later clinical annotation/application, not the opening thesis.

## Known Deviations
- Previous rUTI-led title and early clinical motivation slides are replaced with VF-first framing.
- Main Figures 2-4 move out of the spoken story except where used as appendices.
- AMR is restricted to one appendix backup slide.
`;

await fs.writeFile(path.join(WORKSPACE, "template-audit.txt"), audit, "utf8");

const deviation = `# Deviation Log

## Intentional Scientific Reframe
- Slides 1-10 now present virulence-factor repertoire, curation, and longitudinal stability/change across repeated urinary E. coli isolates.
- Slides 11-12 introduce UTI versus Not_UTI as exploratory clinical annotation/application.
- Appendix slides preserve current UTI-transition and AMR context as backup evidence only.

## Media And Diagram Policy
- Explanatory visuals are rebuilt as native shapes/text.
- Analytical figures are embedded from current PNG outputs.
- No SVG, .bin media, or external image relationships are permitted in the final package.

## Guardrails
- No legacy ASB-vs-UTI conclusion is presented.
- No significant lpf claim is presented.
- No global AMR association claim is presented.
- Stable VF profiles are framed as supportive, not mechanistic proof.
`;

await fs.writeFile(path.join(WORKSPACE, "deviation-log.txt"), deviation, "utf8");

console.log(OUT);

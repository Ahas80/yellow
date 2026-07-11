import fs from "node:fs/promises";

const tmp = "/var/folders/fp/wwwk1rbj70l0k92s92kd2z500000gp/T/codex-presentations/manual-20260710-lecturer-methods/tmp";
const inspectPath = `${tmp}/template-inspect/template-inspect.ndjson`;
const lines = (await fs.readFile(inspectPath, "utf8"))
  .trim()
  .split(/\n+/)
  .map((line) => JSON.parse(line));

const plan = [
  [1, 1, "Minimal claim-linked analysis framing"],
  [2, 3, "Glossary and unit ladder"],
  [3, 4, "Clinical UTI rule"],
  [4, 15, "Two-track clinical and assembly methods flowchart"],
  [5, 5, "Mixed executed set versus Longcycler sensitivity"],
  [6, 6, "Parallel genomic analysis branches"],
  [7, 5, "Longitudinal reconstruction and fresh SNP result"],
  [8, 13, "Supported conclusions and lecturer decisions"],
  [9, 3, "Appendix denominator ladder"],
  [10, 14, "Appendix tool and threshold table"],
  [11, 6, "Appendix distance-method distinctions"],
  [12, 13, "Appendix statistical cautions"],
  [13, 15, "Appendix claim-to-script traceability map"],
  [14, 13, "Appendix assumptions, missing provenance and excluded claims"],
];

const sourceSlides = [...new Set(plan.map((row) => row[1]))];
const sourceLayouts = new Map();
for (const sourceSlide of sourceSlides) {
  const file = `${tmp}/template-inspect/layouts/source-slide-${String(sourceSlide).padStart(2, "0")}.layout.json`;
  sourceLayouts.set(sourceSlide, JSON.parse(await fs.readFile(file, "utf8")));
}

const outputSlides = [];
for (const [outputSlide, sourceSlide, narrativeRole] of plan) {
  const textIds = lines
    .filter((item) => item.slide === sourceSlide && item.kind === "textbox")
    .map((item) => item.id);
  const repositionOrders = new Set(
    sourceSlide === 15
      ? [8, 11, 18, 25, 31, 38, 45, 51, 58, 65]
      : sourceSlide === 1
        ? []
        : [8],
  );
  const repositionIds = sourceLayouts
    .get(sourceSlide)
    .elements
    .filter((element) => typeof element.text === "string" && repositionOrders.has(element.order))
    .map((element) => element.aid);
  const rewriteIds = textIds.filter((id) => !repositionIds.includes(id));
  outputSlides.push({
    outputSlide,
    sourceSlide,
    narrativeRole,
    reuseMode: "duplicate-slide",
    editTargets: [
      {
        action: "rewrite",
        sourceElementIds: rewriteIds,
        reason: "Replace inherited audience copy while preserving the reference layout and styling",
      },
      ...(repositionIds.length
        ? [{
            action: "rewrite-and-reposition",
            sourceElementIds: repositionIds,
            reason: "Shift inherited title/footer text clear of source badges or slide edge after copy replacement",
          }]
        : []),
    ],
  });
}

const allSlides = lines.filter((item) => item.kind === "slide").map((item) => item.slide);
const omittedSourceSlides = allSlides
  .filter((slide) => !sourceSlides.includes(slide))
  .map((sourceSlide) => ({ sourceSlide, reason: "Layout pattern not required for this methods narrative" }));

await fs.writeFile(
  `${tmp}/template-frame-map.json`,
  JSON.stringify({ outputSlides, omittedSourceSlides }, null, 2),
);

const audit = [
  "Template audit",
  "Source: Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx",
  "Purpose: preserve its teaching-first hierarchy, typography, colour system, rails and footers.",
  "Base presentation.pptx was reviewed as a content inventory only because the user explicitly said it did not need to be imported.",
  "Reusable source patterns: question-led opener (1), unit ladder (3), clinical rule (4), repeated episode comparison (5), parallel branches (6), three-card decision close (13), operational map (14), numbered pipeline map (15).",
  "All inherited text boxes on mapped slides are rewrite targets. Slide 4 rewrites the inherited flow-card components in place and replaces only the connector logic to create an accurate two-track flowchart.",
].join("\n");
await fs.writeFile(`${tmp}/template-audit.txt`, audit);

const deviations = [
  "Deviation log",
  "The Base presentation is not cloned; its assembly-to-analysis sequence informs content only.",
  "The compact onboarding deck is the sole visual reference because it is the established readable project style.",
  "Card titles on slide 13 move right within inherited cards; footer text on slides 2-14 moves right to prevent edge clipping.",
  "Slide 4 rewrites the nine inherited flow cards in place and replaces the source arrows because the original connector sequence implied scientifically incorrect single-track attrition. The inherited positions, card components, typography, palette, rails and footer are retained; only connector logic and explanatory labels change.",
  "Slides 6 and 11 remove inherited arrows where methods are parallel rather than sequential.",
  "No external images or additional visual systems are introduced.",
].join("\n");
await fs.writeFile(`${tmp}/deviation-log.txt`, deviations);

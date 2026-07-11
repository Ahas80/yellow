#!/usr/bin/env node

import { createRequire } from "module";
import fs from "fs";
import fsp from "fs/promises";
import path from "path";
import { spawnSync } from "child_process";

const require = createRequire(import.meta.url);
const pptxgen = require("/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/pptxgenjs");
const sharp = require("/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp");

const WORKSPACE = "/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v5";
const OUTPUT = path.join(WORKSPACE, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx");
const PREVIEW_DIR = path.join(WORKSPACE, "preview", "editable");
const QA_DIR = path.join(WORKSPACE, "qa");
const ASSET_DIR = path.join(WORKSPACE, "assets");
const CONTACT_SHEET = path.join(WORKSPACE, "preview", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_contact_sheet.png");
const SKILL_DIR = "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.521.10419/skills/presentations";
const PYTHON = "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";

const ASSETS = {
  vfTopGenes: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_top_gene_prevalence.png",
  moduleGeneCounts: "/Users/Aamir/Desktop/rUTIs/plots/vf/module_gene_counts.png",
  withinHostJaccard: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_within_host_jaccard_distribution.png",
  geneGainLoss: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_gene_gain_loss_consecutive_pairs.png",
  jaccardSameSt: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_jaccard_same_vs_different_st.png",
  geneModelEvidence: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_gene_screening_vs_model_evidence.png",
  geneModelEvidenceSlideCrop: path.join(ASSET_DIR, "vf_gene_screening_vs_model_evidence_slide11_crop.png"),
  main2: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Main_Figure_2_not_uti_to_uti_mechanism_casebook.png",
  main4: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Main_Figure_4_global_vf_signal_and_robustness.png",
  pcoaSt: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_pcoa_jaccard_ST.png",
};

const C = {
  ink: "111827",
  muted: "64748B",
  line: "CBD5E1",
  pale: "F8FAFC",
  white: "FFFFFF",
  blue: "2B6CB0",
  blueFill: "E8F1FA",
  orange: "D97706",
  orangeFill: "FFF4E6",
  green: "2F855A",
  greenFill: "EAF7EF",
  rust: "B65A3C",
  rustFill: "FBEDE8",
  slateFill: "EEF3F8",
  red: "B91C1C",
};

const W = 13.333333;
const H = 7.5;
let pptx;
let SH;

function noLine() {
  return { color: "FFFFFF", transparency: 100 };
}

function solidLine(color = C.line, width = 1) {
  return { color, width };
}

function rect(slide, x, y, w, h, fill = C.white, line = noLine(), radius = false) {
  slide.addShape(radius ? SH.roundRect : SH.rect, {
    x, y, w, h,
    fill: { color: fill },
    line,
    radius: radius ? 0.12 : undefined,
  });
}

function line(slide, x, y, w, h, color = C.line, width = 1) {
  slide.addShape(SH.line, { x, y, w, h, line: { color, width } });
}

function arrow(slide, x1, y1, x2, y2, color = "94A3B8", width = 1.6) {
  slide.addShape(SH.line, {
    x: x1,
    y: y1,
    w: x2 - x1,
    h: y2 - y1,
    line: { color, width, endArrowType: "triangle" },
  });
}

function text(slide, value, x, y, w, h, opts = {}) {
  slide.addText(String(value), {
    x, y, w, h,
    fontFace: opts.face || "Aptos",
    fontSize: opts.size || 14,
    bold: Boolean(opts.bold),
    italic: Boolean(opts.italic),
    color: opts.color || C.ink,
    align: opts.align || "left",
    valign: opts.valign || "top",
    margin: opts.margin ?? 0.04,
    fit: opts.fit || "shrink",
    breakLine: opts.breakLine,
    paraSpaceAfterPt: opts.paraAfter ?? 0,
    breakLine: opts.breakLine,
  });
}

function base(slide, n, eyebrow, title, subtitle = "", opts = {}) {
  slide.background = { color: C.white };
  rect(slide, 0, 0, W, H, C.white);
  rect(slide, 0, 0, 0.13, H, opts.rail || C.blue);
  text(slide, eyebrow.toUpperCase(), 0.56, 0.30, 9.2, 0.22, { size: 8.5, bold: true, color: C.muted, face: "Aptos" });
  text(slide, title, 0.56, 0.58, 11.9, opts.titleH || 0.48, { size: opts.titleSize || 22, bold: true, color: C.ink, face: "Aptos Display" });
  if (subtitle) text(slide, subtitle, 0.56, opts.subtitleY || 1.08, 11.7, 0.32, { size: 11.5, color: C.muted });
  line(slide, 0.56, 1.47, 12.18, 0, C.line, 1);
  line(slide, 0.56, 7.04, 12.18, 0, C.line, 0.8);
  text(slide, "Current pipeline outputs | VF analysis of longitudinal urinary E. coli isolates", 0.56, 7.12, 8.9, 0.18, { size: 7.6, color: C.muted });
  text(slide, String(n).padStart(2, "0"), 11.95, 7.12, 0.78, 0.18, { size: 7.6, color: C.muted, align: "right" });
}

function card(slide, x, y, w, h, title, body, opts = {}) {
  rect(slide, x, y, w, h, opts.fill || C.pale, solidLine(opts.line || C.line, opts.lineW || 0.8), true);
  if (opts.kicker) text(slide, opts.kicker, x + 0.18, y + 0.14, w - 0.36, 0.18, { size: opts.kickerSize || 7.6, bold: true, color: opts.kickerColor || C.muted });
  const ty = y + (opts.kicker ? 0.38 : 0.16);
  text(slide, title, x + 0.18, ty, w - 0.36, opts.titleH || 0.36, { size: opts.titleSize || 13, bold: true, color: opts.titleColor || C.ink });
  if (body) text(slide, body, x + 0.18, opts.bodyY || (ty + 0.38), w - 0.36, h - (opts.bodyY ? (opts.bodyY - y) : 0.58), { size: opts.bodySize || 10.1, color: opts.bodyColor || C.muted, paraAfter: 0 });
}

function metric(slide, x, y, w, value, label, fill, color = C.ink) {
  rect(slide, x, y, w, 0.86, fill, noLine(), true);
  text(slide, value, x + 0.08, y + 0.14, w - 0.16, 0.30, { size: 19, bold: true, color, align: "center" });
  text(slide, label, x + 0.10, y + 0.52, w - 0.20, 0.22, { size: 8.5, color, align: "center" });
}

function label(slide, x, y, w, value, fill = C.slateFill, color = C.muted) {
  rect(slide, x, y, w, 0.28, fill, noLine(), true);
  text(slide, value, x + 0.06, y + 0.06, w - 0.12, 0.14, { size: 7.5, bold: true, color, align: "center" });
}

async function imageContain(slide, imagePath, x, y, w, h, opts = {}) {
  const meta = await sharp(imagePath).metadata();
  const iw = meta.width || 1600;
  const ih = meta.height || 900;
  const ar = iw / ih;
  let dw = w;
  let dh = dw / ar;
  if (dh > h) {
    dh = h;
    dw = dh * ar;
  }
  const dx = x + (w - dw) / 2;
  const dy = y + (h - dh) / 2;
  if (opts.frame !== false) rect(slide, x - 0.04, y - 0.04, w + 0.08, h + 0.08, C.white, solidLine("E5E7EB", 0.8));
  slide.addImage({ path: imagePath, x: dx, y: dy, w: dw, h: dh, altText: opts.alt || path.basename(imagePath) });
}

async function imageCover(slide, imagePath, x, y, w, h, opts = {}) {
  const meta = await sharp(imagePath).metadata();
  const iw = meta.width || 1600;
  const ih = meta.height || 900;
  const ar = iw / ih;
  let dw = w;
  let dh = dw / ar;
  if (dh < h) {
    dh = h;
    dw = dh * ar;
  }
  const dx = x + (w - dw) / 2;
  const dy = y + (h - dh) / 2;
  slide.addImage({ path: imagePath, x: dx, y: dy, w: dw, h: dh, altText: opts.alt || path.basename(imagePath) });
}

function miniLegend(slide, x, y) {
  const items = [
    ["0072B2", "Model not FDR-significant"],
    ["CC79A7", "Model sparse/separation risk"],
    ["D55E00", "Nominal Fisher p < 0.05"],
  ];
  items.forEach((item, i) => {
    const dx = x + i * 2.12;
    line(slide, dx, y + 0.09, 0.30, 0, item[0], 1.4);
    slide.addShape(SH.ellipse, {
      x: dx + 0.12,
      y: y + 0.03,
      w: 0.12,
      h: 0.12,
      fill: { color: item[0] },
      line: noLine(),
    });
    text(slide, item[1], dx + 0.36, y + 0.01, 1.70, 0.16, { size: 7.1, color: C.muted });
  });
}

function slideTitle() {
  const slide = pptx.addSlide();
  slide.background = { color: C.white };
  rect(slide, 0, 0, W, H, C.white);
  rect(slide, 0, 0, 0.16, H, C.blue);
  text(slide, "CURRENT PIPELINE REVIEW / LONGITUDINAL VF ANALYSIS", 0.62, 0.58, 9.4, 0.22, { size: 9, bold: true, color: C.muted });
  text(slide, "Virulence-factor profiling of longitudinal urinary E. coli isolates", 0.62, 1.00, 10.6, 0.86, { size: 31, bold: true, color: C.ink, face: "Aptos Display" });
  text(slide, "Compact onboarding review: current pipeline, core longitudinal VF results, and cautious clinical-status overlay", 0.64, 1.90, 10.2, 0.34, { size: 13, color: C.muted });
  card(slide, 0.64, 2.96, 3.36, 1.36, "What VF repertoire is present?", "Describe the virulence-factor gene space and curated modules in 556 VF-ready urinary E. coli isolates.", { fill: C.blueFill, line: "B9D5F0", kicker: "Question 1", kickerColor: C.blue, bodySize: 9.4 });
  card(slide, 4.40, 2.96, 3.36, 1.36, "How stable are profiles over time?", "Use consecutive within-resident comparisons to quantify VF similarity and gene gain/loss.", { fill: C.greenFill, line: "BCE4C9", kicker: "Question 2", kickerColor: C.green, bodySize: 9.4 });
  card(slide, 8.16, 2.96, 3.36, 1.36, "How should clinical status be overlaid?", "Treat UTI versus Not_UTI as a sparse exploratory annotation, not the whole study frame.", { fill: C.orangeFill, line: "F4C98A", kicker: "Question 3", kickerColor: C.orange, bodySize: 9.4 });
  metric(slide, 1.82, 5.12, 1.75, "556", "VF-ready episodes", C.blueFill, C.blue);
  metric(slide, 4.02, 5.12, 1.75, "227", "VF gene columns", C.slateFill, C.ink);
  metric(slide, 6.22, 5.12, 1.75, "394", "longitudinal pairs", C.greenFill, C.green);
  metric(slide, 8.42, 5.12, 1.75, "17", "UTI annotations", C.orangeFill, C.orange);
  text(slide, "Designed for a 20-22 min presentation plus discussion", 0.62, 7.12, 8.0, 0.18, { size: 8, color: C.muted });
}

function slideYellowLoop() {
  const slide = pptx.addSlide();
  base(slide, 2, "Project orientation", "The YELLOW routine in one plain-English loop", "One episode becomes one bacterial genome profile, then repeated profiles are compared over time.");
  const nodes = [
    ["Clinical episode", "Urine-sampling moment with clinical context.", C.orangeFill, C.orange],
    ["Urinary E. coli isolate", "The bacterium selected for sequencing.", C.blueFill, C.blue],
    ["Genome assembly", "WGS reads become a genome representation.", C.blueFill, C.blue],
    ["VF gene row", "Detected or not detected across VF genes.", C.greenFill, C.green],
    ["Profile outputs", "Modules, burden, similarity, gain/loss.", C.greenFill, C.green],
    ["Clinical overlay", "Status labels added cautiously.", C.orangeFill, C.orange],
  ];
  const x0 = 0.54;
  const boxW = 1.78;
  const gap = 0.28;
  nodes.forEach((node, i) => {
    const x = x0 + i * (boxW + gap);
    card(slide, x, 2.22, boxW, 1.54, node[0], node[1], { fill: node[2], line: node[3], titleColor: C.ink, titleSize: 11.2, titleH: 0.46, bodySize: 8.5, bodyY: 3.00 });
    text(slide, String(i + 1), x + 0.14, 2.36, 0.25, 0.20, { size: 9, bold: true, color: node[3] });
    if (i < nodes.length - 1) arrow(slide, x + boxW + 0.03, 2.98, x + boxW + gap - 0.05, 2.98);
  });
  card(slide, 0.92, 4.82, 3.05, 0.92, "What we measure", "VF gene presence/absence in each selected urinary E. coli assembly.", { fill: C.blueFill, line: "B9D5F0", titleColor: C.blue, titleSize: 11.5, bodySize: 8.8 });
  card(slide, 5.08, 4.82, 3.05, 0.92, "What we compare", "Repeated VF profiles within residents: similarity, stability, and gain/loss.", { fill: C.greenFill, line: "BCE4C9", titleColor: C.green, titleSize: 11.5, bodySize: 8.8 });
  card(slide, 9.24, 4.82, 3.05, 0.92, "Where status enters", "UTI versus Not_UTI is a later exploratory clinical annotation.", { fill: C.orangeFill, line: "F4C98A", titleColor: C.orange, titleSize: 11.5, bodySize: 8.8 });
}

function slideUnitLadder() {
  const slide = pptx.addSlide();
  base(slide, 3, "Project orientation", "Unit changes explain denominator changes", "Each count answers a different question: episode, assembly/profile, feature, pair, or transition.");
  card(slide, 0.76, 1.80, 11.82, 0.76, "Plain rule for reading the deck", "Before interpreting a number, ask what the row represents. Counts change when the analysis moves from clinical records to selected genomes, VF features, repeated pairs, and transition cases.", { fill: C.white, line: C.line, titleSize: 12.4, bodySize: 9.5, bodyY: 2.15 });
  const rows = [
    ["Clinical episodes", "583", "visits/events with primary clinical status", C.orangeFill, C.orange],
    ["VF-ready episodes", "556", "selected episode-level assemblies with VF evidence", C.blueFill, C.blue],
    ["VF feature space", "227 -> 32", "gene columns grouped into curated modules", C.greenFill, C.green],
    ["Longitudinal pairs", "394", "consecutive within-resident comparisons", C.greenFill, C.green],
    ["Clinical transitions", "11", "Not_UTI -> UTI application cases", C.slateFill, C.ink],
  ];
  rows.forEach((row, i) => {
    const y = 3.02 + i * 0.62;
    card(slide, 1.06, y, 2.60, 0.42, row[0], "", { fill: row[3], line: row[4], titleColor: row[4], titleSize: 10.8, titleH: 0.20 });
    text(slide, row[1], 4.22, y - 0.02, 1.60, 0.36, { size: 18, bold: true, color: row[4] });
    text(slide, row[2], 5.88, y + 0.05, 5.40, 0.20, { size: 9.7, color: C.muted });
  });
  label(slide, 8.28, 6.24, 2.90, "Same pipeline, different analysis units", C.slateFill, C.muted);
}

function slideUtiDefinition() {
  const slide = pptx.addSlide();
  base(slide, 4, "Clinical annotation", "How an episode becomes UTI in this project", "Primary UTI status requires culture support plus catheter-aware compatible symptoms.");
  card(slide, 0.86, 1.90, 3.16, 1.56, "1. Culture support", "Urine culture evidence supports possible infection under the primary lower-threshold rule.", { fill: C.blueFill, line: "B9D5F0", titleColor: C.blue, bodySize: 9.7 });
  text(slide, ">=10^3 CFU/mL where CFU data are available", 1.04, 3.04, 2.75, 0.20, { size: 8.2, color: C.blue, bold: true, align: "center" });
  arrow(slide, 4.22, 2.70, 5.06, 2.70);
  text(slide, "+", 4.55, 2.55, 0.28, 0.20, { size: 17, bold: true, color: C.ink, align: "center" });
  card(slide, 5.28, 1.90, 3.16, 1.56, "2. Symptom rule", "Compatible symptoms are interpreted with different logic for catheter and non-catheter episodes.", { fill: C.greenFill, line: "BCE4C9", titleColor: C.green, bodySize: 9.7 });
  text(slide, "catheter-aware symptom definition", 5.70, 3.04, 2.30, 0.20, { size: 8.2, color: C.green, bold: true, align: "center" });
  arrow(slide, 8.62, 2.70, 9.46, 2.70);
  card(slide, 9.68, 2.12, 2.44, 1.10, "UTI", "Both criteria met", { fill: C.orangeFill, line: "F4B15F", titleColor: C.orange, titleSize: 19, bodySize: 9.5 });
  card(slide, 1.02, 4.36, 5.18, 1.04, "Primary clinical denominator", "583 clinical episodes = 18 UTI + 565 Not_UTI", { fill: C.orangeFill, line: "F4B15F", titleColor: C.orange, titleSize: 12.5, bodySize: 12.3, bodyColor: C.ink });
  card(slide, 7.00, 4.36, 5.18, 1.04, "Not_UTI is a mixed comparator", "Episodes that do not meet both criteria; includes heterogeneous clinical and culture contexts.", { fill: C.slateFill, line: C.line, titleColor: C.ink, titleSize: 12.5, bodySize: 9.5 });
  label(slide, 4.12, 6.05, 1.55, "Not_UTI", C.blueFill, C.blue);
  arrow(slide, 5.76, 6.19, 6.62, 6.19, C.line);
  label(slide, 6.74, 6.05, 1.55, "UTI", C.orangeFill, C.orange);
}

function slideDataAsset() {
  const slide = pptx.addSlide();
  base(slide, 5, "Data asset", "VF-ready data link repeated urinary E. coli genomes to clinical episodes", "Dataset-level counts are VF-ready denominators; no disputed clinical participant total is displayed.");
  metric(slide, 0.84, 1.92, 1.76, "556", "episode-level selected assemblies", C.blueFill, C.blue);
  metric(slide, 3.04, 1.92, 1.76, "162", "participants in VF-ready dataset", C.slateFill, C.ink);
  metric(slide, 5.24, 1.92, 1.76, "227", "binary VF gene columns", C.greenFill, C.green);
  metric(slide, 7.44, 1.92, 1.76, "32", "curated modules", C.orangeFill, C.orange);
  metric(slide, 9.64, 1.92, 1.76, "18", "UPEC-candidate modules", C.rustFill, C.rust);
  text(slide, "Repeated-isolate design", 1.02, 3.48, 3.00, 0.24, { size: 15, bold: true, color: C.ink });
  const xs = [1.54, 3.74, 5.94, 8.14, 10.34];
  xs.forEach((x, i) => {
    rect(slide, x, 4.28, 0.66, 0.66, i === 3 ? C.orangeFill : C.blueFill, solidLine(i === 3 ? "F4B15F" : "9AC1E6", 1.4), true);
    text(slide, i === 3 ? "UTI" : "Not_UTI", x - 0.10, 4.50, 0.86, 0.16, { size: 7.4, bold: true, color: i === 3 ? C.orange : C.blue, align: "center" });
    text(slide, `Episode ${i + 1}`, x - 0.18, 5.14, 1.02, 0.16, { size: 8.2, color: C.muted, align: "center" });
    if (i < xs.length - 1) arrow(slide, x + 0.78, 4.61, xs[i + 1] - 0.14, 4.61);
  });
  card(slide, 1.02, 6.02, 11.10, 0.66, "Interpretation stance", "The unit of VF measurement is an episode-level isolate profile. Status labels can be overlaid, but repeated measures, lineage structure, and sparse UTI counts limit causal status claims.", { fill: C.pale, line: C.line, titleSize: 9.7, bodySize: 8.5, bodyY: 6.30 });
}

function slideFeatureFramework() {
  const slide = pptx.addSlide();
  base(slide, 6, "VF representation", "VF features: genes, curated modules, and longitudinal similarity", "Module and score outputs help interpretation; they are not validated disease-causality scores.");
  card(slide, 0.82, 1.92, 2.90, 1.34, "Binary VF gene matrix", "227 presence/absence gene columns per episode-level isolate.", { fill: C.blueFill, line: "B9D5F0", titleColor: C.blue, bodySize: 9.5 });
  card(slide, 5.22, 1.92, 2.90, 1.34, "Curated modules", "32 biological curation units, including 18 UPEC-candidate modules.", { fill: C.greenFill, line: "BCE4C9", titleColor: C.green, bodySize: 9.5 });
  card(slide, 9.62, 1.92, 2.90, 1.34, "Outputs", "Burden, prevalence, Jaccard similarity, gain/loss, and exploratory models.", { fill: C.orangeFill, line: "F4C98A", titleColor: C.orange, bodySize: 9.5 });
  arrow(slide, 3.90, 2.58, 5.02, 2.58);
  arrow(slide, 8.30, 2.58, 9.42, 2.58);
  const matrixX = 1.06;
  const matrixY = 4.18;
  for (let r = 0; r < 6; r += 1) {
    for (let c = 0; c < 10; c += 1) {
      rect(slide, matrixX + c * 0.19, matrixY + r * 0.17, 0.14, 0.12, (r + c * 2) % 3 !== 0 ? C.blue : "E2E8F0");
    }
  }
  text(slide, "Rows are episode-level isolates; columns are VF genes detected from current selected assemblies.", 0.82, 5.42, 3.00, 0.38, { size: 8.8, color: C.muted });
  ["Adhesion", "Iron acquisition", "Secretion", "Toxin", "Capsule", "Unassigned"].forEach((name, i) => {
    const x = 5.26 + (i % 2) * 1.50;
    const y = 4.04 + Math.floor(i / 2) * 0.48;
    label(slide, x, y, 1.26, name, i === 5 ? C.slateFill : C.greenFill, i === 5 ? C.muted : C.green);
  });
  card(slide, 9.38, 3.96, 3.02, 1.34, "Interpretation caution", "Unassigned genes are 25.1% of the VF matrix. Interpret total burden separately from curated and UPEC-candidate summaries.", { fill: "FFF7ED", line: "FED7AA", titleColor: C.orange, bodySize: 9.3 });
  label(slide, 9.58, 5.70, 2.04, "Descriptive framework", C.slateFill, C.muted);
}

async function slideRepertoireModules() {
  const slide = pptx.addSlide();
  base(slide, 7, "VF repertoire", "VF repertoire is summarised as genes and curated modules", "Prevalence and module counts describe what is present; they are not disease-association claims.");
  await imageContain(slide, ASSETS.vfTopGenes, 0.68, 1.74, 5.80, 4.80, { alt: "Top VF gene prevalence" });
  await imageContain(slide, ASSETS.moduleGeneCounts, 6.72, 1.74, 5.82, 4.80, { alt: "VF module gene counts" });
  label(slide, 0.92, 6.55, 1.22, "Descriptive", C.blueFill, C.blue);
  label(slide, 6.94, 6.55, 1.76, "Descriptive framework", C.greenFill, C.green);
}

async function figureSlide(n, eyebrow, title, subtitle, image, callout, evidence, imgBox = {}) {
  const slide = pptx.addSlide();
  base(slide, n, eyebrow, title, subtitle, { titleSize: imgBox.titleSize || 20.5 });
  await imageContain(slide, image, imgBox.x || 0.66, imgBox.y || 1.66, imgBox.w || 9.55, imgBox.h || 4.92, { alt: title });
  if (callout) {
    card(slide, callout.x, callout.y, callout.w, callout.h, callout.title, callout.body, callout);
  }
  if (evidence) {
    label(slide, evidence.x || 10.55, evidence.y || 1.74, evidence.w || 1.80, evidence.text, evidence.fill || C.slateFill, evidence.color || C.muted);
  }
}

async function slideWithinHost() {
  await figureSlide(
    8,
    "Longitudinal stability",
    "Most repeated within-resident VF profiles are highly stable",
    "394 consecutive comparisons from 144 participants; median Jaccard similarity 1.000; 62.4% no VF change.",
    ASSETS.withinHostJaccard,
    {
      x: 10.42,
      y: 1.94,
      w: 1.88,
      h: 1.86,
      title: "Central finding",
      body: "The dominant pattern is VF profile conservation across repeated urinary E. coli isolates from the same resident.",
      fill: C.greenFill,
      line: "BCE4C9",
      titleColor: C.green,
      bodySize: 8.8,
    },
    { text: "Descriptive longitudinal", fill: C.greenFill, color: C.green, w: 2.05 },
    { w: 9.45, h: 4.88 },
  );
}

async function slideGainLoss() {
  await figureSlide(
    9,
    "Longitudinal change",
    "Gain/loss summaries flag candidates for follow-up",
    "Observed VF changes may reflect replacement, assembly/calling variation, or true gene-content change.",
    ASSETS.geneGainLoss,
    {
      x: 10.84,
      y: 1.94,
      w: 1.40,
      h: 1.86,
      title: "Read with context",
      body: "Use lineage and genome-distance evidence before interpreting gain/loss.",
      fill: C.slateFill,
      line: C.line,
      titleColor: C.ink,
      titleSize: 10.5,
      bodySize: 8.2,
    },
    { x: 10.84, y: 4.18, text: "Descriptive", fill: C.slateFill, color: C.muted, w: 1.40 },
    { x: 0.60, y: 1.66, w: 10.10, h: 4.94 },
  );
}

async function slideLineage() {
  await figureSlide(
    10,
    "Lineage context",
    "Sequence-type consistency helps interpret VF stability",
    "Same-ST comparisons support persistent-lineage context, but ST agreement alone does not prove same strain.",
    ASSETS.jaccardSameSt,
    {
      x: 10.42,
      y: 1.94,
      w: 1.88,
      h: 1.92,
      title: "Diagnostic use",
      body: "Lineage context is essential before interpreting profile change or clinical-status contrast.",
      fill: C.slateFill,
      line: C.line,
      titleColor: C.ink,
      bodySize: 8.8,
    },
    { text: "Diagnostic / descriptive", fill: C.slateFill, color: C.muted, w: 2.05 },
    { w: 9.45, h: 4.88 },
  );
}

async function slideClinicalOverlay() {
  const slide = pptx.addSlide();
  base(
    slide,
    11,
    "Clinical annotation",
    "Clinical-status VF signals remain exploratory",
    "The status overlay is useful for hypothesis generation, but it is not a confirmed global VF association result.",
    { titleSize: 20.5 },
  );
  label(slide, 0.72, 1.68, 2.28, "Exploratory phenotype overlay", C.orangeFill, C.orange);
  metric(slide, 3.36, 1.52, 1.52, "17", "UTI episodes", C.orangeFill, C.orange);
  metric(slide, 5.10, 1.52, 1.52, "539", "Not_UTI episodes", C.blueFill, C.blue);
  metric(slide, 6.84, 1.52, 1.52, "0", "FDR-significant global VF hits", C.slateFill, C.ink);
  metric(slide, 8.58, 1.52, 1.52, "37", "sparse/separation flags", C.rustFill, C.rust);

  card(
    slide,
    0.82,
    3.10,
    3.18,
    1.48,
    "1. Initial screen",
    "Some VF genes can look different in simple UTI versus Not_UTI screens.",
    { fill: C.orangeFill, line: "F4C98A", titleColor: C.orange, titleSize: 12.5, bodySize: 9.6 },
  );
  label(slide, 1.15, 4.16, 1.16, "nominal only", C.white, C.orange);
  card(
    slide,
    5.08,
    3.10,
    3.18,
    1.48,
    "2. Model check",
    "After participant-aware modelling and FDR correction, no global VF association remains significant.",
    { fill: C.slateFill, line: C.line, titleColor: C.ink, titleSize: 12.5, bodySize: 9.4 },
  );
  label(slide, 5.44, 4.16, 1.46, "0 corrected hits", C.white, C.ink);
  card(
    slide,
    9.34,
    3.10,
    3.18,
    1.48,
    "3. Interpretation",
    "The sparse UTI denominator and repeated/lineage structure make this a prioritisation layer, not a discovery claim.",
    { fill: C.greenFill, line: "BCE4C9", titleColor: C.green, titleSize: 12.5, bodySize: 9.2 },
  );
  label(slide, 9.70, 4.16, 1.76, "hypothesis-generating", C.white, C.green);

  arrow(slide, 4.16, 3.84, 4.88, 3.84, "94A3B8", 1.5);
  arrow(slide, 8.42, 3.84, 9.14, 3.84, "94A3B8", 1.5);

  rect(slide, 1.08, 5.42, 11.00, 0.82, C.white, solidLine("E5E7EB", 0.8), true);
  text(slide, "Take-home message", 1.36, 5.68, 1.72, 0.18, { size: 10.2, bold: true, color: C.ink });
  text(slide, "The clinical overlay helps decide which VF features or cases deserve follow-up, but the current population-level evidence does not support a robust UTI-associated VF signature.", 3.10, 5.62, 8.50, 0.28, { size: 10.2, color: C.muted });
}

function slideParticipant20026() {
  const slide = pptx.addSlide();
  base(slide, 12, "Clinical application", "Participant 20026: stable VF profile despite symptom emergence", "This example applies the VF pipeline to a clinical transition; it does not prove mechanism.");
  line(slide, 1.24, 3.20, 10.04, 0, "94A3B8", 2.2);
  card(slide, 1.74, 2.70, 1.62, 0.92, "T3", "Not_UTI", { fill: C.blueFill, line: "9AC1E6", titleColor: C.blue, titleSize: 13, bodySize: 10.5, bodyColor: C.blue });
  card(slide, 9.02, 2.70, 1.62, 0.92, "UTI-1", "UTI", { fill: C.orangeFill, line: "F4B15F", titleColor: C.orange, titleSize: 13, bodySize: 10.5, bodyColor: C.orange });
  arrow(slide, 3.50, 3.20, 8.88, 3.20);
  label(slide, 5.50, 2.80, 1.06, "42 days", C.slateFill, C.ink);
  metric(slide, 2.00, 4.40, 1.58, "5", "SNPs", C.greenFill, C.green);
  metric(slide, 4.08, 4.40, 1.58, "+0/-0", "VF genes", C.greenFill, C.green);
  metric(slide, 6.16, 4.40, 1.58, "stable", "VF modules", C.greenFill, C.green);
  metric(slide, 8.24, 4.40, 1.58, "strong", "same-strain evidence", C.greenFill, C.green);
  card(slide, 0.84, 5.78, 5.52, 0.72, "Transition casebook context", "11 clinical Not_UTI -> UTI transitions; 10 WGS/VF-linked; 1 missing endpoint.", { fill: C.pale, line: C.line, titleSize: 9.5, bodySize: 8.3, bodyY: 6.10 });
  card(slide, 6.78, 5.78, 5.52, 0.72, "Mechanism summary", "4 stable-profile transitions and 3 replacement-consistent transitions; evidence buckets organise cases, not proof.", { fill: C.pale, line: C.line, titleSize: 9.5, bodySize: 8.3, bodyY: 6.10 });
}

function slideTakeaways() {
  const slide = pptx.addSlide();
  base(slide, 13, "Close", "What this VF-first review establishes and what remains uncertain", "Use these as discussion prompts and handover priorities.");
  card(slide, 0.86, 1.92, 3.42, 2.74, "Evidence", "The current pipeline produces a coherent VF-ready longitudinal dataset: 556 isolate profiles, 227 VF gene columns, 32 modules, and 394 consecutive within-resident comparisons.", { fill: C.blueFill, line: "B9D5F0", titleColor: C.blue, bodySize: 10.8 });
  card(slide, 4.96, 1.92, 3.42, 2.74, "Boundary", "VF profiles are often stable, but presence/absence does not measure expression or activity. Clinical status associations remain exploratory.", { fill: C.orangeFill, line: "F4C98A", titleColor: C.orange, bodySize: 10.8 });
  card(slide, 9.06, 1.92, 3.42, 2.74, "Next steps", "Use outputs to prioritise lineage-aware longitudinal follow-up, expression/regulation hypotheses, and targeted review of transition cases.", { fill: C.greenFill, line: "BCE4C9", titleColor: C.green, bodySize: 10.8 });
  text(slide, "Discussion prompt: what additional evidence would separate stable carriage plus host-state change from unmeasured bacterial regulation?", 1.06, 5.58, 11.20, 0.34, { size: 14, bold: true, color: C.ink, align: "center" });
}

function slideHandover() {
  const slide = pptx.addSlide();
  base(slide, 14, "Appendix / handover", "Practical map: where to enter and rerun the VF pipeline", "Keep this as the operational handover slide for a new project member.");
  const left = [
    ["Clinical key", "results/clinical/status_map.csv"],
    ["VF matrix", "results/vf/vf_pa_all.csv"],
    ["Analysis dataset", "results/vf/vf_analysis_ready.csv"],
    ["Modules", "results/vf/gene_module_map.csv"],
    ["Longitudinal", "results/vf/vf_longitudinal_transitions.csv"],
    ["Final figures", "35_final_figure_pack.R"],
  ];
  const right = [
    ["Counts / validation", "results/final_figures/final_figure_validation_checks.csv"],
    ["Denominator flow", "results/qc/pipeline_denominator_summary.csv"],
    ["VF diagnostics", "results/vf/vf_dataset_diagnostics.txt"],
    ["VF visual audit", "results/vf/vf_visualisation_audit.csv"],
    ["Casebook", "results/mechanism/not_uti_to_uti_casebook.md"],
    ["Architecture", "docs/pipeline_architecture.md"],
  ];
  text(slide, "Canonical outputs", 0.86, 1.72, 2.6, 0.24, { size: 14, bold: true, color: C.blue });
  text(slide, "Sources of truth", 6.82, 1.72, 2.6, 0.24, { size: 14, bold: true, color: C.green });
  const render = (rows, x) => rows.forEach((r, i) => {
    const y = 2.08 + i * 0.58;
    rect(slide, x, y, 5.28, 0.46, i % 2 ? C.white : C.pale, solidLine("E5E7EB", 0.6), true);
    text(slide, r[0], x + 0.16, y + 0.10, 1.40, 0.16, { size: 8.4, bold: true, color: C.ink });
    text(slide, r[1], x + 1.62, y + 0.10, 3.42, 0.16, { size: 7.3, color: C.muted, face: "Aptos Mono" });
  });
  render(left, 0.86);
  render(right, 6.82);
  card(slide, 2.02, 6.12, 9.30, 0.54, "Rerun order", "Clinical status -> assembly QC/canonical selection -> VF matrix -> VF analysis dataset -> modules/scores -> longitudinal summaries -> final figure pack.", { fill: C.slateFill, line: C.line, titleSize: 9.2, bodySize: 8.2, bodyY: 6.36 });
}

function slideDetailedPipeline() {
  const slide = pptx.addSlide();
  base(slide, 15, "Appendix / detailed map", "Numbered clinical-to-VF pipeline with validated handoff counts", "This is the dense version for behind-the-scenes questions.");
  const stages = [
    ["1", "Clinical classification", "00a / 00b", "585 classified episodes before manual exclusions"],
    ["2", "Primary inclusion", "status_map.csv", "583 included: 18 UTI, 565 Not_UTI; 2 excluded"],
    ["3", "Assembly QC", "12a_wgs_qc.R", "1,291 QC records incl. assembler alternatives"],
    ["4", "Canonical selection", "canonical_assembly_selection.csv", "556 selected episode-level assemblies"],
    ["5", "VF P/A matrix", "02_gene_presence_analysis.R", "556 episode-level VF rows"],
    ["6", "VF/model-ready", "22_vf_build_analysis_dataset.R", "556 episodes: 17 UTI, 539 Not_UTI; 27 clinical episodes lack VF-ready evidence"],
    ["7", "Feature framework", "26_vf_define_gene_modules.R", "227 VF genes -> 32 modules, incl. 18 UPEC-candidate modules"],
    ["8", "Longitudinal pairs", "24_vf_longitudinal_dynamics.R", "394 consecutive comparisons from 144 participants"],
    ["9", "Clinical application", "casebook outputs", "11 Not_UTI -> UTI transitions; 10 WGS/VF-linked; 1 missing endpoint"],
  ];
  stages.forEach((s, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    const x = 0.74 + col * 4.06;
    const y = 1.78 + row * 1.54;
    const fill = i < 2 ? C.orangeFill : i < 6 ? C.blueFill : i === 7 ? C.greenFill : C.slateFill;
    card(slide, x, y, 3.52, 1.08, s[1], s[3], { fill, line: "D8DEE8", titleSize: 10.4, bodySize: 8.1, bodyY: y + 0.62 });
    rect(slide, x + 0.12, y + 0.12, 0.28, 0.28, C.ink, noLine(), true);
    text(slide, s[0], x + 0.12, y + 0.17, 0.28, 0.10, { size: 6.8, bold: true, color: C.white, align: "center" });
    text(slide, s[2], x + 0.52, y + 0.40, 2.80, 0.14, { size: 6.9, color: C.muted, face: "Aptos Mono" });
    if (col < 2) arrow(slide, x + 3.56, y + 0.54, x + 3.96, y + 0.54, "94A3B8", 1.2);
  });
}

async function slideCasebook() {
  await figureSlide(
    16,
    "Appendix / transition evidence",
    "Full transition mechanism casebook",
    "Use for questions about Not_UTI -> UTI transitions; buckets organise evidence, not proof.",
    ASSETS.main2,
    null,
    { x: 10.18, y: 1.74, text: "Descriptive casebook", fill: C.slateFill, color: C.muted, w: 2.02 },
    { x: 0.64, y: 1.66, w: 11.20, h: 4.94, titleSize: 21 },
  );
}

async function slideRobustnessLineage() {
  const slide = pptx.addSlide();
  base(slide, 17, "Appendix / interpretation checks", "Robustness and lineage diagnostics stay as backup interpretation aids", "Use these for questions about exploratory clinical signals and lineage structure.");
  await imageContain(slide, ASSETS.main4, 0.66, 1.74, 5.78, 4.74, { alt: "Main Figure 4" });
  await imageContain(slide, ASSETS.pcoaSt, 6.72, 1.74, 5.78, 4.74, { alt: "VF Jaccard PCoA by ST" });
  label(slide, 0.92, 6.55, 2.32, "No global VF FDR signal", C.orangeFill, C.orange);
  label(slide, 7.00, 6.55, 2.20, "Lineage diagnostic", C.slateFill, C.muted);
}

async function build() {
  for (const [name, p] of Object.entries(ASSETS)) {
    if (name !== "geneModelEvidenceSlideCrop" && !fs.existsSync(p)) throw new Error(`Missing asset ${name}: ${p}`);
  }
  await fsp.mkdir(path.dirname(OUTPUT), { recursive: true });
  await fsp.mkdir(PREVIEW_DIR, { recursive: true });
  await fsp.mkdir(QA_DIR, { recursive: true });
  await fsp.mkdir(ASSET_DIR, { recursive: true });

  await sharp(ASSETS.geneModelEvidence)
    .extract({ left: 0, top: 150, width: 3300, height: 1500 })
    .png()
    .toFile(ASSETS.geneModelEvidenceSlideCrop);

  pptx = new pptxgen();
  SH = pptx.ShapeType || {};
  SH.rect = SH.rect || "rect";
  SH.roundRect = SH.roundRect || "roundRect";
  SH.line = SH.line || "line";
  pptx.defineLayout({ name: "CUSTOM_WIDE", width: W, height: H });
  pptx.layout = "CUSTOM_WIDE";
  pptx.author = "Aamir / Codex";
  pptx.subject = "Compact editable VF onboarding deck";
  pptx.title = "Longitudinal urinary E. coli VF pipeline review compact onboarding";
  pptx.company = "rUTI project";
  pptx.lang = "en-US";
  pptx.theme = {
    headFontFace: "Aptos Display",
    bodyFontFace: "Aptos",
    lang: "en-US",
  };

  slideTitle();
  slideYellowLoop();
  slideUnitLadder();
  slideUtiDefinition();
  slideDataAsset();
  slideFeatureFramework();
  await slideRepertoireModules();
  await slideWithinHost();
  await slideGainLoss();
  await slideLineage();
  await slideClinicalOverlay();
  slideParticipant20026();
  slideTakeaways();
  slideHandover();
  slideDetailedPipeline();
  await slideCasebook();
  await slideRobustnessLineage();

  await pptx.writeFile({ fileName: OUTPUT });

  const qa = spawnSync(PYTHON, [path.join(WORKSPACE, "scripts", "qa_pptx_package.py"), OUTPUT], { encoding: "utf8" });
  await fsp.writeFile(path.join(QA_DIR, "pptx_package_check.txt"), `${qa.stdout || ""}${qa.stderr || ""}`);
  if (qa.status !== 0) {
    throw new Error(`PPTX package QA failed:\n${qa.stdout}\n${qa.stderr}`);
  }

  console.log(OUTPUT);
}

build().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

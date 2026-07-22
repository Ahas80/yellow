#!/usr/bin/env node

import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";

import {
  createSlideContext,
  ensureArtifactToolWorkspace,
  importArtifactTool,
  padSlideNumber,
  saveBlobToFile,
} from "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.709.11516/skills/presentations/container_tools/artifact_tool_utils.mjs";

import {
  DEFAULT_PROJECT_ROOT,
  PYTHON,
  SKILL_DIR,
  assertNoEmptyPlaceholders,
  assertPptxHasNoForbiddenContent,
  copyFileAtomic,
  loadRegistry,
  parseArgs,
  releaseCounts,
  run,
  writeJsonAtomic,
  writeTextAtomic,
} from "../../../../../scripts/generated/longcycler_release_presentation_common.mjs";

const SIZE = { width: 1280, height: 720 };
const MIN_BODY_FONT_SIZE = 16;
const CANONICAL_ROOT = path.join(
  DEFAULT_PROJECT_ROOT,
  "outputs",
  "manual-20260527-current-review",
  "presentations",
);

const C = {
  ink: "#111827",
  muted: "#64748B",
  light: "#F8FAFC",
  line: "#CBD5E1",
  blue: "#2B6CB0",
  blueFill: "#E8F1FA",
  orange: "#D97706",
  orangeFill: "#FFF4E6",
  green: "#2F855A",
  greenFill: "#EAF7EF",
  rust: "#B65A3C",
  rustFill: "#FBEDE8",
  slateFill: "#EEF3F8",
  purple: "#7C3AED",
  cyan: "#0891B2",
  white: "#FFFFFF",
};

const VARIANTS = {
  v3: {
    slides: 22,
    family: "ruti-longitudinal-vf-pipeline-review-v3",
    output: "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx",
    footer: "Current Longcycler-only pipeline outputs | operational UTI phenotype",
    audience: "Scientific review",
  },
  v4: {
    slides: 25,
    family: "ruti-longitudinal-vf-pipeline-review-v4",
    output: "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx",
    footer: "Longcycler-only review with onboarding | operational UTI phenotype",
    audience: "Scientific review with onboarding",
  },
  "v5-full": {
    slides: 17,
    family: "ruti-longitudinal-vf-pipeline-review-v5",
    output: "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-28.pptx",
    footer: "Longcycler-only current review | operational UTI phenotype",
    audience: "Full 17-slide review",
  },
  "v5-compact": {
    slides: 17,
    family: "ruti-longitudinal-vf-pipeline-review-v5",
    output: "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx",
    footer: "Longcycler-only compact onboarding | operational UTI phenotype",
    audience: "Compact onboarding review",
  },
};

const PLOT_CANDIDATES = {
  topGenes: ["vf_top_gene_prevalence.png"],
  modules: ["module_gene_counts.png"],
  stability: ["vf_within_host_jaccard_distribution.png"],
  gainLoss: ["vf_gene_gain_loss_consecutive_pairs.png"],
  lineageBox: ["vf_jaccard_same_vs_different_st.png"],
  clinical: ["vf_gene_screening_vs_model_evidence.png"],
  moduleStatus: ["module_prevalence_by_status.png"],
  casebook: ["not_uti_to_uti_case_matrix.png"],
  stableState: ["host_context_transition_heatmap.png"],
  robustness: ["model_stability_flags.png"],
  lineage: ["vf_pcoa_jaccard_ST.png"],
  accessory: ["FigS07_plasmid_amr_context.png"],
};

function requireVariant(value) {
  if (!VARIANTS[value]) {
    throw new Error(`--variant must be one of ${Object.keys(VARIANTS).join(", ")}`);
  }
  return value;
}

async function walk(root) {
  const out = [];
  for (const entry of await fs.readdir(root, { withFileTypes: true }).catch(() => [])) {
    const full = path.join(root, entry.name);
    if (entry.isDirectory()) out.push(...(await walk(full)));
    else if (entry.isFile()) out.push(full);
  }
  return out;
}

async function resolvePlots(registry, projectRoot) {
  const registered = (registry.plot_files || []).map((item) => path.resolve(String(item)));
  const inventory = registered.length
    ? registered
    : await walk(path.join(projectRoot, "plots"));
  const available = inventory.filter((item) => /\.(?:png|jpe?g)$/i.test(item));
  const resolved = {};
  for (const [key, candidates] of Object.entries(PLOT_CANDIDATES)) {
    const match = available.find((item) => candidates.includes(path.basename(item)));
    if (match) resolved[key] = match;
  }
  const required = ["topGenes", "modules", "stability", "gainLoss", "lineageBox", "clinical", "moduleStatus", "casebook", "stableState", "robustness", "lineage"];
  const missing = required.filter((key) => !resolved[key]);
  if (missing.length) {
    throw new Error(`Missing current-run plot assets for: ${missing.join(", ")}`);
  }
  return resolved;
}

function setNotes(slide, notes) {
  if (!notes) return;
  if (slide.speakerNotes?.textFrame?.setText) {
    slide.speakerNotes.textFrame.setText(notes);
    slide.speakerNotes.setVisible?.(true);
    return;
  }
  if (slide.speakerNotes?.setText) slide.speakerNotes.setText(notes);
}

function shape(ctx, slide, frame, options = {}) {
  return ctx.addShape(slide, {
    ...frame,
    geometry: options.geometry || "rect",
    fill: options.fill || C.white,
    line: options.line || ctx.line(),
    name: options.name,
  });
}

function text(ctx, slide, frame, value, options = {}) {
  const role = options.role || "body";
  const fontSize = Number(options.size ?? MIN_BODY_FONT_SIZE);
  if (role === "body" && fontSize < MIN_BODY_FONT_SIZE) {
    throw new Error(`Body copy must be at least ${MIN_BODY_FONT_SIZE}pt: ${String(value).slice(0, 80)}`);
  }
  return ctx.addText(slide, {
    ...frame,
    text: String(value),
    fontSize,
    bold: options.bold || false,
    color: options.color || C.ink,
    align: options.align || "left",
    valign: options.valign || "top",
    typeface: options.typeface || (options.mono ? "Aptos Mono" : "Aptos"),
    insets: options.insets || { left: 0, right: 0, top: 0, bottom: 0 },
    name: options.name || `${role}-${Math.round(frame.x)}-${Math.round(frame.y)}-${Math.round(frame.w)}-${Math.round(frame.h)}`,
  });
}

function base(ctx, slide, spec, variant, slideNumber) {
  shape(ctx, slide, { x: 0, y: 0, w: ctx.W, h: ctx.H }, { fill: C.white, name: "background" });
  shape(ctx, slide, { x: 0, y: 0, w: 12, h: ctx.H }, { fill: spec.rail || C.blue, name: "accent-rail" });
  text(ctx, slide, { x: 54, y: 32, w: 1040, h: 18 }, (spec.eyebrow || "CURRENT REVIEW").toUpperCase(), { size: 11.5, bold: true, color: C.muted, name: "eyebrow", role: "chrome" });
  text(ctx, slide, { x: 54, y: 57, w: 1160, h: 44 }, spec.title, { size: spec.titleSize || 29, bold: true, typeface: "Aptos Display", name: "title", role: "heading" });
  if (spec.subtitle) text(ctx, slide, { x: 54, y: 105, w: 1140, h: 28 }, spec.subtitle, { size: MIN_BODY_FONT_SIZE, color: C.muted, name: "subtitle" });
  shape(ctx, slide, { x: 54, y: 140, w: 1172, h: 1.2 }, { fill: C.line, name: "header-rule" });
  shape(ctx, slide, { x: 54, y: 674, w: 1172, h: 1.1 }, { fill: C.line, name: "footer-rule" });
  text(ctx, slide, { x: 54, y: 682, w: 1040, h: 18 }, VARIANTS[variant].footer, { size: 9.5, color: C.muted, name: "footer", role: "chrome" });
  text(ctx, slide, { x: 1150, y: 682, w: 76, h: 18 }, String(slideNumber).padStart(2, "0"), { size: 9.5, color: C.muted, align: "right", name: "slide-number", role: "chrome" });
}

function addCard(ctx, slide, frame, titleValue, bodyValue, options = {}) {
  shape(ctx, slide, frame, {
    fill: options.fill || C.light,
    line: { style: "solid", fill: options.line || C.line, width: 1 },
    geometry: "roundRect",
    name: options.name || `card-${titleValue}`,
  });
  let top = frame.y + 16;
  if (options.kicker) {
    text(ctx, slide, { x: frame.x + 16, y: top, w: frame.w - 32, h: 18 }, options.kicker, { size: 12, bold: true, color: options.kickerColor || C.muted, role: "label" });
    top += 24;
  }
  text(ctx, slide, { x: frame.x + 16, y: top, w: frame.w - 32, h: options.titleH || 26 }, titleValue, { size: options.titleSize || 17, bold: true, color: options.titleColor || C.ink, role: "heading" });
  text(ctx, slide, { x: frame.x + 16, y: top + (options.titleH || 26) + 6, w: frame.w - 32, h: frame.y + frame.h - top - (options.titleH || 26) - 20 }, bodyValue, { size: Math.max(MIN_BODY_FONT_SIZE, Number(options.bodySize ?? MIN_BODY_FONT_SIZE)), color: options.bodyColor || C.muted });
}

function addMetric(ctx, slide, frame, value, labelValue, options = {}) {
  shape(ctx, slide, frame, { fill: options.fill || C.slateFill, geometry: "roundRect", name: `metric-${labelValue}` });
  text(ctx, slide, { x: frame.x + 10, y: frame.y + 12, w: frame.w - 20, h: 31 }, value, { size: options.valueSize || 27, bold: true, color: options.color || C.ink, align: "center", role: "heading" });
  text(ctx, slide, { x: frame.x + 10, y: frame.y + 50, w: frame.w - 20, h: frame.h - 54 }, labelValue, { size: Math.max(MIN_BODY_FONT_SIZE, Number(options.labelSize ?? MIN_BODY_FONT_SIZE)), color: options.color || C.ink, align: "center" });
}

function addLabel(ctx, slide, frame, value, options = {}) {
  shape(ctx, slide, frame, { fill: options.fill || C.slateFill, geometry: "roundRect", name: `label-${value}` });
  text(ctx, slide, { x: frame.x + 8, y: frame.y + 5, w: frame.w - 16, h: frame.h - 8 }, value, { size: options.size || 12, bold: true, color: options.color || C.ink, align: "center", role: "label" });
}

function addArrow(ctx, slide, from, to, color = "#94A3B8", name = "connector") {
  const dx = to.x - from.x;
  const dy = to.y - from.y;
  const length = Math.max(1, Math.sqrt(dx * dx + dy * dy));
  const angle = Math.atan2(dy, dx) * 180 / Math.PI;
  const line = shape(ctx, slide, { x: (from.x + to.x) / 2 - length / 2, y: (from.y + to.y) / 2 - 1, w: length, h: 2 }, { fill: color, name });
  line.position.rotation = angle;
  const head = shape(ctx, slide, { x: to.x - 7, y: to.y - 5, w: 10, h: 10 }, { fill: color, geometry: "triangle", name: `${name}-head` });
  head.position.rotation = angle + 90;
}

async function addFigure(ctx, slide, imagePath, frame, alt) {
  shape(ctx, slide, { x: frame.x - 7, y: frame.y - 7, w: frame.w + 14, h: frame.h + 14 }, { fill: C.white, line: { style: "solid", fill: "#E5E7EB", width: 1 }, name: `figure-frame-${path.basename(imagePath)}` });
  await ctx.addImage(slide, { ...frame, path: imagePath, fit: "contain", alt, name: `figure-${path.basename(imagePath)}` });
}

function titleSlide(ctx, slide, variant, c) {
  shape(ctx, slide, { x: 0, y: 0, w: ctx.W, h: ctx.H }, { fill: C.white, name: "background" });
  shape(ctx, slide, { x: 0, y: 0, w: 16, h: ctx.H }, { fill: C.blue, name: "accent-rail" });
  text(ctx, slide, { x: 62, y: 54, w: 1080, h: 22 }, "CURRENT REVIEW  /  SELECTED LONGCYCLER ANALYSIS", { size: 12, bold: true, color: C.muted, role: "chrome" });
  text(ctx, slide, { x: 62, y: 92, w: 1120, h: 90 }, "Virulence-factor profiling of longitudinal urinary E. coli isolates", { size: 39, bold: true, typeface: "Aptos Display" });
  text(ctx, slide, { x: 62, y: 190, w: 1000, h: 28 }, "One selected QC-passing assembly per episode; operational UTI overlay remains exploratory", { size: 17, color: C.muted });
  const cards = [
    ["Question 1", "What VF repertoire is present?", `${c.episodes} selected episodes; ${c.vfFeatures} binary VF features.`, C.blueFill, C.blue],
    ["Question 2", "How stable are profiles over time?", `${c.directPairs} direct pairs; ${c.transitions} adjacent transitions.`, C.greenFill, C.green],
    ["Question 3", "How should clinical status be overlaid?", `${c.uti} operational UTI and ${c.notUti} operational Not_UTI episodes.`, C.orangeFill, C.orange],
  ];
  cards.forEach((item, index) => addCard(ctx, slide, { x: 64 + index * 398, y: 278, w: 340, h: 145 }, item[1], item[2], { kicker: item[0], fill: item[3], line: item[3], kickerColor: item[4], titleColor: item[4], bodySize: 12.5 }));
  const metrics = [
    [String(c.episodes), "selected episodes", C.blueFill, C.blue],
    [String(c.residents), "residents", C.slateFill, C.ink],
    [String(c.transitions), "adjacent transitions", C.greenFill, C.green],
    [String(c.uti), "operational UTI", C.orangeFill, C.orange],
  ];
  metrics.forEach((item, index) => addMetric(ctx, slide, { x: 190 + index * 225, y: 500, w: 185, h: 86 }, item[0], item[1], { fill: item[2], color: item[3] }));
  text(ctx, slide, { x: 62, y: 682, w: 1020, h: 18 }, VARIANTS[variant].footer, { size: 9.5, color: C.muted, role: "chrome" });
  text(ctx, slide, { x: 1150, y: 682, w: 76, h: 18 }, "01", { size: 9.5, color: C.muted, align: "right", role: "chrome" });
}

function scopeSlide(ctx, slide, variant, number, c, registry) {
  base(ctx, slide, {
    eyebrow: "Scope lock",
    title: "The analytical scope is selected QC-passing Longcycler only",
    subtitle: "Clinical status is an operational annotation; the scientific object is the repeated isolate profile.",
  }, variant, number);
  const y = 190;
  addCard(ctx, slide, { x: 76, y, w: 330, h: 300 }, "Selected episode-level genomes", `${c.episodes} episodes from ${c.residents} residents contribute one selected QC-passing assembly per episode. The selected manifest is the only analytical genome source.`, { fill: C.blueFill, line: "#BDD7EE", titleColor: C.blue, bodySize: 14.5 });
  addCard(ctx, slide, { x: 475, y, w: 330, h: 300 }, "Repeated-measures design", `${c.directPairs} direct within-resident pairs are available for broad profile comparisons. ${c.transitions} adjacent transitions from ${c.transitionResidents} residents define the temporal transition layer.`, { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 14.5 });
  addCard(ctx, slide, { x: 874, y, w: 330, h: 300 }, "Operational clinical overlay", `${c.uti} episodes are operational UTI and ${c.notUti} are operational Not_UTI. Interpret status comparisons as ${registry.analysis_scope.interpretation}.`, { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 14.5 });
  addLabel(ctx, slide, { x: 418, y: 536, w: 442, h: 31 }, "No alternate assembly source enters the analytical cohort", { fill: C.slateFill, color: C.ink, size: 11.5 });
}

function loopSlide(ctx, slide, variant, number) {
  base(ctx, slide, {
    eyebrow: "Project orientation",
    title: "The review follows one plain-English loop",
    subtitle: "One episode becomes one selected genome profile, then repeated profiles are compared over time.",
  }, variant, number);
  const steps = [
    ["Clinical episode", "Apply the operational clinical definition", C.orangeFill, C.orange],
    ["Selected genome", "Use the QC-passing Longcycler assembly", C.blueFill, C.blue],
    ["VF gene row", "Record binary VF feature presence", C.greenFill, C.green],
    ["Longitudinal pair", "Compare profiles within resident", C.slateFill, C.ink],
    ["Clinical overlay", "Interpret status cautiously", C.orangeFill, C.orange],
  ];
  for (let i = 0; i < steps.length - 1; i += 1) addArrow(ctx, slide, { x: 248 + i * 236, y: 292 }, { x: 288 + i * 236, y: 292 }, "#94A3B8", `loop-${i + 1}`);
  steps.forEach((item, index) => addCard(ctx, slide, { x: 54 + index * 236, y: 224, w: 194, h: 137 }, item[0], item[1], { fill: item[2], line: item[2], titleColor: item[3], bodySize: 11.5, titleSize: 15 }));
  addCard(ctx, slide, { x: 112, y: 438, w: 300, h: 105 }, "What we measure", "VF gene presence/absence in each selected episode-level genome.", { fill: C.blueFill, line: "#BDD7EE", titleSize: 15, bodySize: 12.5 });
  addCard(ctx, slide, { x: 490, y: 438, w: 300, h: 105 }, "What we compare", "Repeated profiles within residents, directly and as adjacent transitions.", { fill: C.greenFill, line: "#BCE4C9", titleSize: 15, bodySize: 12.5 });
  addCard(ctx, slide, { x: 868, y: 438, w: 300, h: 105 }, "Where status enters", "Only after the genomic and temporal denominators are explicit.", { fill: C.orangeFill, line: "#F4C98A", titleSize: 15, bodySize: 12.5 });
}

function unitsSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Project orientation",
    title: "Unit changes explain denominator changes",
    subtitle: "Each count answers a different question: episode, resident, profile, pair, or transition.",
  }, variant, number);
  const rows = [
    ["Source clinical episodes", String(c.sourceEpisodes), `${c.sourceResidents} residents; retained for attrition/QC only`, C.orangeFill, C.orange],
    ["Selected analytical episodes", String(c.episodes), `${c.residents} residents; one selected genome per episode`, C.blueFill, C.blue],
    ["Direct within-resident pairs", String(c.directPairs), "all available within-resident comparisons", C.slateFill, C.ink],
    ["Adjacent transitions", String(c.transitions), `${c.transitionResidents} residents; ${c.transitionsAtThreshold} at or below ${c.snpThreshold} SNPs`, C.greenFill, C.green],
    ["Focused Not_UTI→UTI", String(c.toUti), `${c.toUtiAtThreshold} at or below ${c.snpThreshold} SNPs`, C.orangeFill, C.orange],
  ];
  rows.forEach((item, index) => {
    const y = 175 + index * 85;
    addLabel(ctx, slide, { x: 150, y, w: 260, h: 34 }, item[0], { fill: item[3], color: item[4], size: 11.5 });
    text(ctx, slide, { x: 454, y: y - 2, w: 110, h: 38 }, item[1], { size: 24, bold: true, color: item[4], align: "right" });
    text(ctx, slide, { x: 616, y: y + 2, w: 530, h: 34 }, item[2], { size: MIN_BODY_FONT_SIZE, color: C.muted });
  });
  addLabel(ctx, slide, { x: 448, y: 610, w: 384, h: 31 }, "Always name the unit beside the count", { fill: C.slateFill, color: C.ink, size: 11 });
}

function phenotypeSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Clinical annotation",
    title: "How an episode becomes operational UTI in this project",
    subtitle: "The operational label is a study annotation, not a universal clinical diagnosis.",
  }, variant, number);
  addArrow(ctx, slide, { x: 390, y: 285 }, { x: 456, y: 285 }, "#94A3B8", "phenotype-1");
  addArrow(ctx, slide, { x: 794, y: 285 }, { x: 860, y: 285 }, "#94A3B8", "phenotype-2");
  addCard(ctx, slide, { x: 82, y: 210, w: 308, h: 150 }, "1. Culture support", "The episode passes the project’s culture and isolate evidence rules.", { fill: C.blueFill, line: "#BDD7EE", titleColor: C.blue, bodySize: 13.5 });
  addCard(ctx, slide, { x: 456, y: 210, w: 338, h: 150 }, "2. Symptom rule", "Compatible symptoms are interpreted with the operational definition and available episode context.", { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 13.5 });
  addCard(ctx, slide, { x: 860, y: 230, w: 270, h: 110 }, "Operational UTI", `${c.uti} selected analytical episodes`, { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 13.5 });
  addCard(ctx, slide, { x: 106, y: 440, w: 490, h: 105 }, "Selected analytical denominator", `${c.episodes} episodes: ${c.uti} operational UTI and ${c.notUti} operational Not_UTI.`, { fill: C.orangeFill, line: "#F4C98A", titleSize: 15, bodySize: 13 });
  addCard(ctx, slide, { x: 682, y: 440, w: 490, h: 105 }, "Not_UTI is a mixed comparator", "It means the episode did not meet the operational UTI label; it is not a single biological state.", { fill: C.slateFill, line: C.line, titleSize: 15, bodySize: 13 });
}

function cohortSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Data asset",
    title: "VF-ready data link selected genomes to repeated clinical episodes",
    subtitle: "All analytical denominators are selected Longcycler outputs; source counts appear only as QC context.",
  }, variant, number);
  const metrics = [
    [c.episodes, "selected episodes", C.blueFill, C.blue],
    [c.residents, "residents", C.slateFill, C.ink],
    [c.vfFeatures, "binary VF features", C.greenFill, C.green],
    [c.uti, "operational UTI", C.orangeFill, C.orange],
    [c.notUti, "operational Not_UTI", C.rustFill, C.rust],
  ];
  metrics.forEach((item, index) => addMetric(ctx, slide, { x: 66 + index * 236, y: 176, w: 190, h: 86 }, String(item[0]), item[1], { fill: item[2], color: item[3], valueSize: 25 }));
  text(ctx, slide, { x: 98, y: 326, w: 360, h: 24 }, "Repeated-isolate design", { size: 19, bold: true });
  const xs = [150, 370, 590, 810, 1030];
  for (let i = 0; i < xs.length - 1; i += 1) addArrow(ctx, slide, { x: xs[i] + 80, y: 437 }, { x: xs[i + 1] - 10, y: 437 }, "#94A3B8", `episode-${i + 1}`);
  xs.forEach((x, index) => {
    shape(ctx, slide, { x, y: 400, w: 70, h: 70 }, { fill: index === 3 ? C.orangeFill : C.blueFill, line: { style: "solid", fill: index === 3 ? "#F4B15F" : "#9AC1E6", width: 2 }, geometry: "ellipse", name: `episode-node-${index + 1}` });
    text(ctx, slide, { x: x - 8, y: 426, w: 86, h: 18 }, index === 3 ? "UTI" : "Not_UTI", { size: 10.5, bold: true, color: index === 3 ? C.orange : C.blue, align: "center", role: "label" });
    text(ctx, slide, { x: x - 16, y: 484, w: 102, h: 18 }, `Episode ${index + 1}`, { size: 10.5, color: C.muted, align: "center", role: "label" });
  });
  addCard(ctx, slide, { x: 98, y: 548, w: 1050, h: 72 }, "Interpretation stance", "Episode profiles are repeated within residents. Status comparisons must account for repeated measures, lineage structure, and the sparse operational UTI denominator.", { fill: C.light, titleSize: 13.5, bodySize: 11.5 });
}

function frameworkSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "VF representation",
    title: "VF features are summarised as genes, curated modules, and longitudinal similarity",
    subtitle: "Module and score outputs support interpretation; they are not validated disease-causality measures.",
  }, variant, number);
  addArrow(ctx, slide, { x: 355, y: 257 }, { x: 480, y: 257 }, "#94A3B8", "framework-1");
  addArrow(ctx, slide, { x: 795, y: 257 }, { x: 920, y: 257 }, "#94A3B8", "framework-2");
  addCard(ctx, slide, { x: 74, y: 188, w: 280, h: 138 }, "Binary VF gene matrix", `${c.vfFeatures} presence/absence features per selected episode-level genome.`, { fill: C.blueFill, line: "#B9D5F0", titleColor: C.blue });
  addCard(ctx, slide, { x: 480, y: 188, w: 315, h: 138 }, "Curated biological modules", "Genes are grouped to make repertoire patterns navigable without changing the underlying binary evidence.", { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green });
  addCard(ctx, slide, { x: 920, y: 188, w: 285, h: 138 }, "Longitudinal outputs", "Prevalence, Jaccard similarity, gain/loss summaries, and exploratory clinical overlays.", { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange });
  for (let r = 0; r < 6; r += 1) for (let col = 0; col < 11; col += 1) shape(ctx, slide, { x: 102 + col * 19, y: 392 + r * 17, w: 14, h: 12 }, { fill: (r + col * 2) % 3 ? C.blue : "#E2E8F0", name: `matrix-${r}-${col}` });
  ["Adhesion", "Iron acquisition", "Secretion", "Toxin", "Capsule", "Unassigned"].forEach((name, index) => addLabel(ctx, slide, { x: 500 + (index % 2) * 145, y: 380 + Math.floor(index / 2) * 45, w: 126, h: 29 }, name, { fill: index === 5 ? C.slateFill : C.greenFill, color: index === 5 ? C.muted : C.green, size: 9.5 }));
  addCard(ctx, slide, { x: 900, y: 372, w: 305, h: 145 }, "Interpretation caution", "Gene presence does not establish expression, activity, disease mechanism, or a causal clinical-state effect.", { fill: C.orangeFill, line: "#FED7AA", titleColor: C.orange, bodySize: 13.5 });
  addLabel(ctx, slide, { x: 520, y: 558, w: 240, h: 30 }, "Descriptive feature framework", { fill: C.slateFill, color: C.muted });
}

async function figureSlide(ctx, slide, variant, number, spec) {
  base(ctx, slide, spec, variant, number);
  await addFigure(ctx, slide, spec.image, spec.frame || { x: 62, y: 160, w: 930, h: 472 }, spec.alt || spec.title);
  if (spec.callout) addCard(ctx, slide, spec.callout.frame, spec.callout.title, spec.callout.body, spec.callout.options || {});
  if (spec.label) addLabel(ctx, slide, spec.label.frame, spec.label.text, spec.label.options || {});
}

function pairSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Longitudinal denominator",
    title: "Direct pairs and adjacent transitions answer different questions",
    subtitle: "The transition layer is temporally ordered; the direct-pair layer describes the wider repeated-isolate network.",
  }, variant, number);
  addMetric(ctx, slide, { x: 114, y: 190, w: 240, h: 104 }, String(c.directPairs), "all direct within-resident pairs", { fill: C.blueFill, color: C.blue, valueSize: 31, labelSize: 12 });
  addMetric(ctx, slide, { x: 520, y: 190, w: 240, h: 104 }, String(c.transitions), "adjacent transitions", { fill: C.greenFill, color: C.green, valueSize: 31, labelSize: 12 });
  addMetric(ctx, slide, { x: 926, y: 190, w: 240, h: 104 }, String(c.transitionResidents), "residents with adjacent transitions", { fill: C.slateFill, color: C.ink, valueSize: 31, labelSize: 12 });
  addCard(ctx, slide, { x: 114, y: 370, w: 350, h: 150 }, "Broad profile question", "Use all direct pairs to describe how similar episode-level VF profiles are across the available within-resident history.", { fill: C.blueFill, line: "#BDD7EE", titleColor: C.blue, bodySize: 14 });
  addCard(ctx, slide, { x: 490, y: 370, w: 350, h: 150 }, "Temporal transition question", "Use adjacent transitions when asking what changed between successive observed episodes.", { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 14 });
  addCard(ctx, slide, { x: 866, y: 370, w: 300, h: 150 }, "Genome-distance context", `${c.transitionsAtThreshold} of ${c.transitions} adjacent transitions are at or below the operational ${c.snpThreshold}-SNP threshold.`, { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 14 });
}

function transitionSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Clinical transition application",
    title: "Focused transition evidence is aggregate, linked, and deliberately cautious",
    subtitle: "The release reports the full focused casebook; it does not elevate one anecdote into a mechanism claim.",
  }, variant, number);
  addArrow(ctx, slide, { x: 376, y: 310 }, { x: 442, y: 310 }, "#94A3B8", "transition-1");
  addArrow(ctx, slide, { x: 804, y: 310 }, { x: 870, y: 310 }, "#94A3B8", "transition-2");
  addCard(ctx, slide, { x: 88, y: 226, w: 288, h: 170 }, "Focused Not_UTI→UTI transitions", `${c.toUti} adjacent transitions enter the focused clinical-transition view.`, { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, titleH: 48, bodySize: 14 });
  addCard(ctx, slide, { x: 442, y: 226, w: 362, h: 170 }, `At or below ${c.snpThreshold} SNPs`, `${c.toUtiAtThreshold} of ${c.toUti} focused transitions fall at or below the operational genome-distance threshold.`, { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 14 });
  addCard(ctx, slide, { x: 870, y: 226, w: 310, h: 170 }, "Mechanism casebook", `${c.casebookCases} cases; ${c.casebookLinked} linked; ${c.casebookMissing} missing endpoints.`, { fill: C.blueFill, line: "#BDD7EE", titleColor: C.blue, bodySize: 14 });
  addCard(ctx, slide, { x: 170, y: 492, w: 940, h: 82 }, "Interpretation boundary", "Low genome distance and stable VF content can support a stable-lineage interpretation, but they cannot by themselves identify host-state change, bacterial regulation, or clinical causation.", { fill: C.light, line: C.line, titleSize: 14, bodySize: 12.5 });
}

function closeSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Close",
    title: "What the Longcycler-only review establishes and what remains uncertain",
    subtitle: "Use these as the three discussion anchors.",
  }, variant, number);
  addCard(ctx, slide, { x: 82, y: 184, w: 330, h: 300 }, "Evidence", `${c.episodes} selected profiles from ${c.residents} residents support ${c.directPairs} direct comparisons and ${c.transitions} adjacent transitions. The mechanism casebook is complete (${c.casebookLinked}/${c.casebookCases} linked).`, { fill: C.blueFill, line: "#B9D5F0", titleColor: C.blue, bodySize: 14.5 });
  addCard(ctx, slide, { x: 474, y: 184, w: 330, h: 300 }, "Boundary", `Operational UTI is sparse (${c.uti}/${c.episodes}). VF presence/absence does not measure expression, regulation, pathogenic activity, or a causal clinical-state effect.`, { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 14.5 });
  addCard(ctx, slide, { x: 866, y: 184, w: 330, h: 300 }, "Next steps", "Prioritise lineage-aware longitudinal analyses, targeted transition review, and additional evidence that can separate bacterial change from host-state and sampling effects.", { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 14.5 });
  text(ctx, slide, { x: 100, y: 542, w: 1080, h: 32 }, "Discussion: what extra evidence would distinguish stable carriage plus host-state change from unmeasured bacterial regulation?", { size: 17, bold: true, align: "center" });
}

function handoverSlide(ctx, slide, variant, number) {
  base(ctx, slide, {
    eyebrow: "Appendix / handover",
    title: "Practical map: where to enter and rerun the Longcycler-only pipeline",
    subtitle: "The numbered research-question layer runs from RQ01 through RQ10.",
  }, variant, number);
  const rows = [
    ["Selected manifest", "results/assembly/canonical_assembly_selection.csv", "One QC-passing Longcycler assembly per episode"],
    ["Clinical key", "results/clinical/status_map.csv", "Operational UTI annotation and audit fields"],
    ["VF matrix", "results/vf/vf_pa_all.csv", "Episode-level binary VF profiles"],
    ["Direct pairs", "results/strain/pairwise_metrics.csv", "All within-resident comparisons"],
    ["Adjacent transitions", "results/longitudinal/longcycler_transitions.csv", "Temporally ordered comparisons"],
    ["Focused casebook", "results/mechanism/not_uti_to_uti_casebook.csv", "Nine linked transition cases"],
    ["Research questions", "scripts/research_questions/run_all.R", "RQ01–RQ10 release layer"],
  ];
  rows.forEach((item, index) => {
    const y = 158 + index * 68;
    shape(ctx, slide, { x: 76, y, w: 1080, h: 52 }, { fill: index % 2 ? C.white : C.light, line: { style: "solid", fill: "#E5E7EB", width: 1 }, name: `handover-row-${index + 1}` });
    text(ctx, slide, { x: 96, y: y + 8, w: 170, h: 34 }, item[0], { size: 12.5, bold: true, color: C.blue, role: "label" });
    text(ctx, slide, { x: 286, y: y + 8, w: 410, h: 36 }, item[1], { size: MIN_BODY_FONT_SIZE, mono: true });
    text(ctx, slide, { x: 720, y: y + 8, w: 410, h: 36 }, item[2], { size: MIN_BODY_FONT_SIZE, color: C.muted });
  });
}

function pipelineSlide(ctx, slide, variant, number, c) {
  base(ctx, slide, {
    eyebrow: "Appendix / pipeline",
    title: "Numbered clinical-to-VF pipeline with validated handoff counts",
    subtitle: "RQ01–RQ10 consume the same selected analytical cohort and current-run outputs.",
    titleSize: 27,
  }, variant, number);
  const stages = [
    ["1", "Source/QC context", `${c.sourceEpisodes} episodes; ${c.sourceResidents} residents`, C.orangeFill],
    ["2", "Selected cohort", `${c.episodes} episodes; ${c.residents} residents`, C.blueFill],
    ["3", "Operational phenotype", `${c.uti} UTI; ${c.notUti} Not_UTI`, C.orangeFill],
    ["4", "VF feature matrix", `${c.vfFeatures} binary features`, C.blueFill],
    ["5", "MLST context", `${c.mlstTyped} typed episodes; ${c.mlstDistinct} labels`, C.slateFill],
    ["6", "Direct pairs", `${c.directPairs} comparisons`, C.blueFill],
    ["7", "Adjacent transitions", `${c.transitions} from ${c.transitionResidents} residents`, C.greenFill],
    ["8", "Focused transition view", `${c.toUti} transitions; ${c.toUtiAtThreshold} ≤${c.snpThreshold} SNPs`, C.orangeFill],
    ["9", "Mechanism casebook", `${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing}: cases/linked/missing`, C.greenFill],
    ["10", "Research-question release", "RQ01–RQ10; none retired", C.slateFill],
  ];
  for (let i = 0; i < stages.length; i += 1) {
    const col = i % 5;
    const row = Math.floor(i / 5);
    const x = 54 + col * 238;
    const y = 178 + row * 210;
    if (col < 4) addArrow(ctx, slide, { x: x + 214, y: y + 74 }, { x: x + 230, y: y + 74 }, "#94A3B8", `stage-${i + 1}`);
  }
  stages.forEach((item, index) => {
    const col = index % 5;
    const row = Math.floor(index / 5);
    const x = 54 + col * 238;
    const y = 178 + row * 210;
    shape(ctx, slide, { x, y, w: 214, h: 148 }, { fill: item[3], line: { style: "solid", fill: "#D8DEE8", width: 1 }, geometry: "roundRect", name: `pipeline-${item[0]}` });
    shape(ctx, slide, { x: x + 14, y: y + 16, w: 30, h: 30 }, { fill: C.ink, geometry: "ellipse", name: `pipeline-number-${item[0]}` });
    text(ctx, slide, { x: x + 14, y: y + 22, w: 30, h: 16 }, item[0], { size: 11, bold: true, color: C.white, align: "center", role: "label" });
    text(ctx, slide, { x: x + 54, y: y + 15, w: 145, h: 44 }, item[1], { size: 13.5, bold: true, role: "heading" });
    text(ctx, slide, { x: x + 16, y: y + 72, w: 182, h: 64 }, item[2], { size: MIN_BODY_FONT_SIZE, color: C.muted });
  });
}

function sourceSlide(ctx, slide, variant, number, c, registryPath, registrySha) {
  base(ctx, slide, {
    eyebrow: "Appendix / evidence registry",
    title: "Sources of truth and lookup resources for questions",
    subtitle: "Counts and plots are bound to the current Longcycler release claim registry.",
  }, variant, number);
  const rows = [
    ["Release claim registry", registryPath],
    ["Registry SHA-256", registrySha],
    ["Selected cohort", `${c.episodes} episodes / ${c.residents} residents / ${c.uti} UTI / ${c.notUti} Not_UTI`],
    ["Pair layers", `${c.directPairs} direct pairs; ${c.transitions} adjacent transitions from ${c.transitionResidents} residents`],
    ["Focused casebook", `${c.toUti} transitions; ${c.toUtiAtThreshold} ≤${c.snpThreshold} SNPs; ${c.casebookCases}/${c.casebookLinked}/${c.casebookMissing} cases/linked/missing`],
    ["Near-miss audit", `${c.nearMiss} rows; separate from operational UTI cases`],
    ["Research questions", "RQ01–RQ10"],
  ];
  rows.forEach((item, index) => {
    const y = 160 + index * 65;
    shape(ctx, slide, { x: 82, y, w: 1116, h: 49 }, { fill: index % 2 ? C.white : C.light, line: { style: "solid", fill: "#E5E7EB", width: 1 }, name: `source-row-${index + 1}` });
    text(ctx, slide, { x: 102, y: y + 8, w: 210, h: 30 }, item[0], { size: 12.5, bold: true, color: C.blue, role: "label" });
    text(ctx, slide, { x: 336, y: y + 8, w: 836, h: 34 }, item[1], { size: MIN_BODY_FONT_SIZE, color: index < 2 ? C.muted : C.ink, mono: index < 2 });
  });
  addLabel(ctx, slide, { x: 410, y: 626, w: 460, h: 30 }, "583/166/18/565 is attrition/QC context only", { fill: C.orangeFill, color: C.orange, size: 11 });
}

function notes(titleValue, say, boundary) {
  return [`SLIDE: ${titleValue}`, `SAY: ${say}`, `BOUNDARY: ${boundary}`].join("\n\n");
}

async function buildNarrative({ presentation, artifact, workspace, variant, c, registry, registryPath, registrySha, plots }) {
  const spec = VARIANTS[variant];
  const slides = [];
  const add = async (renderer, note) => {
    const slide = presentation.slides.add();
    const number = slides.length + 1;
    const ctx = createSlideContext(artifact, { slideSize: SIZE, slideNumber: number, workspaceDir: workspace, assetDir: path.join(workspace, "assets"), outputDir: path.join(workspace, "output") });
    await renderer(ctx, slide, number);
    setNotes(slide, note);
    slides.push(slide);
  };

  const standardFigures = {
    top: {
      eyebrow: "VF repertoire", title: "Common VF genes describe the observed repertoire, not disease association", subtitle: `Prevalence across ${c.episodes} selected episode-level genomes.`, image: plots.topGenes,
      frame: { x: 62, y: 158, w: 900, h: 478 }, alt: "Top virulence-factor gene prevalence plot",
      callout: { frame: { x: 1000, y: 190, w: 200, h: 182 }, title: "Read as", body: "Descriptive prevalence ranking only. Gene presence does not imply expression, activity, or UTI association.", options: { fill: C.blueFill, line: "#B9D5F0", titleColor: C.blue, bodySize: 12.5 } },
      label: { frame: { x: 1018, y: 406, w: 164, h: 30 }, text: "Descriptive", options: { fill: C.blueFill, color: C.blue } },
    },
    modules: {
      eyebrow: "VF framework", title: "Curated modules organise genes into interpretable biological systems", subtitle: "Modules are navigation units; they are not validated disease-causality scores.", image: plots.modules,
      frame: { x: 62, y: 158, w: 890, h: 478 }, alt: "Virulence-factor module gene-count plot",
      callout: { frame: { x: 988, y: 188, w: 214, h: 204 }, title: "Why it matters", body: "The module framework structures descriptive repertoire and longitudinal summaries while retaining the binary gene evidence.", options: { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 12.5 } },
    },
    stability: {
      eyebrow: "Longitudinal stability", title: "Repeated within-resident VF profiles are assessed on the explicit pair layers", subtitle: `${c.directPairs} direct pairs; ${c.transitions} adjacent transitions from ${c.transitionResidents} residents.`, image: plots.stability,
      frame: { x: 62, y: 158, w: 910, h: 478 }, alt: "Within-resident virulence-factor Jaccard distribution",
      callout: { frame: { x: 1004, y: 188, w: 198, h: 204 }, title: "Central reading", body: "Profile stability is descriptive evidence. Use adjacent transitions and lineage context before interpreting change.", options: { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 12.5 } },
      label: { frame: { x: 1014, y: 430, w: 178, h: 30 }, text: "Descriptive", options: { fill: C.greenFill, color: C.green } },
    },
    gainLoss: {
      eyebrow: "Longitudinal change", title: "Gain/loss summaries flag candidates for follow-up", subtitle: `Observed changes are reported on ${c.transitions} adjacent transitions and require lineage plus genome-distance context.`, image: plots.gainLoss,
      frame: { x: 58, y: 158, w: 1010, h: 478 }, alt: "Virulence-factor gene gain and loss across adjacent transitions",
      callout: { frame: { x: 1092, y: 188, w: 128, h: 220 }, title: "Caution", body: "Change may reflect lineage replacement, calling differences, or true gene-content change.", options: { fill: C.slateFill, line: C.line, titleSize: 14, bodySize: 11 } },
    },
    lineageBox: {
      eyebrow: "Lineage context", title: "Sequence-type context helps interpret VF similarity", subtitle: `${c.mlstTyped} selected episodes have preferred MLST labels spanning ${c.mlstDistinct} distinct labels.`, image: plots.lineageBox,
      frame: { x: 62, y: 158, w: 900, h: 478 }, alt: "Virulence-factor similarity by sequence-type agreement",
      callout: { frame: { x: 1000, y: 188, w: 202, h: 206 }, title: "Diagnostic use", body: "ST agreement is context, not proof of the same strain. Pairwise genome distance remains essential.", options: { fill: C.slateFill, line: C.line, bodySize: 12.5 } },
    },
    clinical: {
      eyebrow: "Clinical annotation", title: "Clinical-status VF signals remain exploratory", subtitle: `${c.uti} operational UTI and ${c.notUti} operational Not_UTI selected episodes; no causal interpretation.`, image: plots.clinical,
      frame: { x: 62, y: 158, w: 900, h: 478 }, alt: "Exploratory virulence-factor clinical-status model evidence",
      callout: { frame: { x: 1000, y: 188, w: 202, h: 218 }, title: "Safe claim", body: "Nominal screens and model diagnostics generate hypotheses; they do not establish a robust status association.", options: { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 12.5 } },
      label: { frame: { x: 1000, y: 438, w: 202, h: 30 }, text: "Exploratory overlay", options: { fill: C.orangeFill, color: C.orange } },
    },
  };

  await add((ctx, slide) => titleSlide(ctx, slide, variant, c), notes("Title", "Set the scope immediately: selected QC-passing Longcycler only.", "Operational UTI is an exploratory annotation."));

  if (variant === "v4") {
    await add((ctx, slide, n) => loopSlide(ctx, slide, variant, n), notes("Plain-English loop", "Walk through episode, selected genome, VF profile, repeated pair, then clinical overlay.", "Do not collapse episode, resident, pair, and transition units."));
    await add((ctx, slide, n) => unitsSlide(ctx, slide, variant, n, c), notes("Denominator units", "Name the unit beside every count.", "Source counts are QC context only."));
    await add((ctx, slide, n) => phenotypeSlide(ctx, slide, variant, n, c), notes("Operational phenotype", "Explain the project’s operational UTI annotation.", "Not_UTI is a mixed comparator, not a single biological state."));
  }

  if (variant === "v5-full" || variant === "v5-compact") {
    await add((ctx, slide, n) => loopSlide(ctx, slide, variant, n), notes("Plain-English loop", "Orient the audience before showing analytical plots.", "Keep the selected assembly and operational phenotype explicit."));
    await add((ctx, slide, n) => unitsSlide(ctx, slide, variant, n, c), notes("Units", "Explain why counts change when the unit changes.", "583 source episodes are attrition/QC context only."));
    await add((ctx, slide, n) => phenotypeSlide(ctx, slide, variant, n, c), notes("Operational phenotype", "Define how the project labels episodes.", "The label is not a universal clinical diagnosis."));
    await add((ctx, slide, n) => cohortSlide(ctx, slide, variant, n, c), notes("Selected cohort", "Anchor the analytical denominator.", "All analytical genome content comes from selected Longcycler outputs."));
    await add((ctx, slide, n) => frameworkSlide(ctx, slide, variant, n, c), notes("VF representation", "Explain genes, modules, and longitudinal summaries.", "Presence/absence is not expression or function."));
    await add(async (ctx, slide, n) => {
      base(ctx, slide, { eyebrow: "VF framework", title: "VF repertoire is summarised as genes and curated modules", subtitle: "Both plots describe the selected analytical cohort; neither plot is a causal clinical model." }, variant, n);
      await addFigure(ctx, slide, plots.topGenes, { x: 58, y: 170, w: 552, h: 420 }, "Top virulence-factor genes");
      await addFigure(ctx, slide, plots.modules, { x: 670, y: 170, w: 552, h: 420 }, "Virulence-factor modules");
      addLabel(ctx, slide, { x: 128, y: 610, w: 280, h: 28 }, `${c.vfFeatures} binary VF features`, { fill: C.blueFill, color: C.blue });
      addLabel(ctx, slide, { x: 806, y: 610, w: 280, h: 28 }, "Curated module framework", { fill: C.greenFill, color: C.green });
    }, notes("Paired repertoire views", "Show gene-level and module-level views together.", "Treat modules as curation units."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.stability), notes("Stability", "Use the explicit direct and adjacent pair denominators.", "Profile stability is descriptive."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.gainLoss), notes("Gain/loss", "Flag candidates for follow-up.", "Change needs lineage and genome-distance context."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.lineageBox), notes("Lineage", "Use MLST as a diagnostic layer.", "ST agreement does not prove the same strain."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.clinical), notes("Clinical overlay", "Keep the sparse operational phenotype denominator visible.", "Exploratory observational analysis only."));
    await add((ctx, slide, n) => transitionSlide(ctx, slide, variant, n, c), notes("Focused transitions", "Report 9 transitions, 5 at threshold, and a complete 9-case casebook.", "Do not infer mechanism from low SNP distance alone."));
    await add((ctx, slide, n) => closeSlide(ctx, slide, variant, n, c), notes("Close", "Summarise evidence, boundary, and next steps.", "Keep causality outside the supported claim set."));
    await add((ctx, slide, n) => handoverSlide(ctx, slide, variant, n), notes("Handover", "Show the operational entry points.", "Research questions run RQ01 through RQ10."));
    await add((ctx, slide, n) => pipelineSlide(ctx, slide, variant, n, c), notes("Pipeline", "Use the numbered release pipeline for behind-the-scenes questions.", "All stages consume the selected analytical cohort."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, { eyebrow: "Appendix / transition evidence", title: "Full focused transition mechanism casebook", subtitle: `${c.casebookCases} cases; ${c.casebookLinked} linked; ${c.casebookMissing} missing. Buckets organise evidence, not proof.`, image: plots.casebook, frame: { x: 62, y: 158, w: 1110, h: 482 }, alt: "Focused Not_UTI to UTI mechanism casebook" }), notes("Casebook", "Use the complete nine-case casebook for transition questions.", "Buckets organise evidence; they do not prove mechanism."));
    await add(async (ctx, slide, n) => {
      base(ctx, slide, { eyebrow: "Appendix / interpretation checks", title: "Robustness and lineage diagnostics stay as backup interpretation aids", subtitle: `The ${c.nearMiss}-row near-miss audit remains separate from operational UTI cases.` }, variant, n);
      await addFigure(ctx, slide, plots.robustness, { x: 54, y: 170, w: 552, h: 418 }, "Global virulence-factor robustness plot");
      await addFigure(ctx, slide, plots.lineage, { x: 670, y: 170, w: 552, h: 418 }, "Virulence-factor lineage PCoA");
      addLabel(ctx, slide, { x: 166, y: 612, w: 300, h: 28 }, `${c.uti} operational UTI; exploratory`, { fill: C.orangeFill, color: C.orange });
      addLabel(ctx, slide, { x: 802, y: 612, w: 280, h: 28 }, "Lineage diagnostic", { fill: C.slateFill, color: C.muted });
    }, notes("Robustness and lineage", "Keep these two diagnostic aids together.", "Near-miss rows are not operational UTI cases."));
  } else {
    await add((ctx, slide, n) => scopeSlide(ctx, slide, variant, n, c, registry), notes("Scope", "State the selected Longcycler-only analytical lock.", "The phenotype overlay is exploratory."));
    await add((ctx, slide, n) => cohortSlide(ctx, slide, variant, n, c), notes("Selected cohort", "Anchor 532 episodes, 161 residents, and the operational phenotype split.", "Source clinical counts are QC context only."));
    await add((ctx, slide, n) => pipelineSlide(ctx, slide, variant, n, c), notes("Pipeline", "Walk through the current release handoffs.", "RQ01–RQ10 use the same selected cohort."));
    await add((ctx, slide, n) => frameworkSlide(ctx, slide, variant, n, c), notes("VF representation", "Explain the descriptive feature framework.", "Presence/absence is not expression."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.top), notes("VF repertoire", "Read the top-gene ranking descriptively.", "Do not infer clinical association."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.modules), notes("VF modules", "Use modules for navigation.", "Modules are not validated disease scores."));
    await add((ctx, slide, n) => pairSlide(ctx, slide, variant, n, c), notes("Pair layers", "Separate 893 direct pairs from 371 adjacent transitions.", "Name the unit beside each count."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.gainLoss), notes("Gain/loss", "Show adjacent-transition change candidates.", "Changes need lineage and QC context."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.lineageBox), notes("Lineage", "Use MLST as contextual evidence.", "ST agreement is not proof of the same strain."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, standardFigures.clinical), notes("Clinical overlay", "Keep 16 operational UTI and 516 Not_UTI visible.", "Exploratory observational analysis only."));
    await add((ctx, slide, n) => transitionSlide(ctx, slide, variant, n, c), notes("Focused transitions", "Report the aggregate focused transition evidence.", "Avoid anecdotal mechanism claims."));
    await add((ctx, slide, n) => closeSlide(ctx, slide, variant, n, c), notes("Close", "Land on evidence, boundary, and next steps.", "Do not cross the causal boundary."));
    await add((ctx, slide, n) => handoverSlide(ctx, slide, variant, n), notes("Handover", "Show the practical entry points.", "The research-question layer is RQ01–RQ10."));
    await add((ctx, slide, n) => unitsSlide(ctx, slide, variant, n, c), notes("Denominator audit", "Explain that 583/166/18/565 is source/QC context only.", "Analytical claims use 532/161/16/516."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, { eyebrow: "Appendix / clinical annotation", title: "Module prevalence by status is an exploratory annotation view", subtitle: `Operational phenotype denominator: ${c.uti} UTI and ${c.notUti} Not_UTI selected episodes.`, image: plots.moduleStatus, frame: { x: 62, y: 158, w: 1040, h: 480 }, alt: "Virulence-factor module prevalence by operational UTI status", callout: { frame: { x: 1120, y: 190, w: 100, h: 170 }, title: "Label", body: "Exploratory clinical annotation.", options: { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, titleSize: 13, bodySize: 10.5 } } }), notes("Module prevalence", "Use only as an exploratory status view.", "Repeated measures and lineage remain interpretive limits."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, { eyebrow: "Appendix / transition evidence", title: "Full focused transition mechanism casebook", subtitle: `${c.casebookCases} cases; ${c.casebookLinked} linked; ${c.casebookMissing} missing. Buckets organise evidence, not proof.`, image: plots.casebook, frame: { x: 62, y: 158, w: 1110, h: 482 }, alt: "Focused Not_UTI to UTI mechanism casebook" }), notes("Casebook", "Use the full nine-case linked casebook.", "Casebook buckets do not prove mechanism."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, { eyebrow: "Appendix / stable-state context", title: "Stable strain and changing clinical state", subtitle: `${c.toUtiAtThreshold} of ${c.toUti} focused transitions are at or below ${c.snpThreshold} SNPs.`, image: plots.stableState, frame: { x: 62, y: 158, w: 1110, h: 482 }, alt: "Stable strain and host-context interpretation plot" }), notes("Stable-state context", "Use low distance and stable VF as supporting evidence.", "Do not claim host-state mechanism from genomic stability alone."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, { eyebrow: "Appendix / robustness", title: "Population-level robustness boundary", subtitle: `${c.uti} operational UTI episodes limit adjusted inference; the ${c.nearMiss}-row near-miss audit is separate.`, image: plots.robustness, frame: { x: 62, y: 158, w: 1110, h: 482 }, alt: "Global virulence-factor robustness plot" }), notes("Robustness", "Keep the sparse operational phenotype denominator explicit.", "Near-miss rows are not operational UTI cases."));
    await add((ctx, slide, n) => figureSlide(ctx, slide, variant, n, { eyebrow: "Appendix / lineage", title: "Lineage structure is an interpretation check", subtitle: "Virulence-factor profile PCoA by preferred sequence-type label.", image: plots.lineage, frame: { x: 68, y: 158, w: 1020, h: 482 }, alt: "Virulence-factor PCoA by sequence type", callout: { frame: { x: 1102, y: 190, w: 118, h: 190 }, title: "Use for", body: "Questions about lineage structure and VF clustering.", options: { fill: C.slateFill, line: C.line, titleSize: 13, bodySize: 10.5 } } }), notes("Lineage PCoA", "Use as a diagnostic check before status interpretation.", "Clustering does not establish causality."));
    await add(async (ctx, slide, n) => {
      base(ctx, slide, { eyebrow: "Appendix / accessory context", title: "Accessory and mobile-element context remains secondary", subtitle: "Use only when it helps interpret a transition; this review does not make a dedicated AMR or causal accessory claim." }, variant, n);
      if (plots.accessory) await addFigure(ctx, slide, plots.accessory, { x: 62, y: 158, w: 930, h: 480 }, "Accessory mobile-element transition context");
      else addCard(ctx, slide, { x: 104, y: 210, w: 740, h: 250 }, "No dedicated accessory figure in the current plot registry", "The release keeps accessory and mobile-element evidence as contextual backup. It does not substitute that evidence for VF, lineage, or pairwise genome-distance interpretation.", { fill: C.slateFill, line: C.line, titleSize: 21, bodySize: 16 });
      addCard(ctx, slide, { x: 1010, y: 190, w: 190, h: 218 }, "Guardrail", "No dedicated AMR association claim. No causal interpretation. Use transition-level context only.", { fill: C.rustFill, line: "#EDB9A7", titleColor: C.rust, bodySize: 12.5 });
    }, notes("Accessory context", "Keep accessory information in a supporting role.", "No dedicated AMR or causal claim."));
    await add((ctx, slide, n) => sourceSlide(ctx, slide, variant, n, c, registryPath, registrySha), notes("Sources", "Use the claim registry for provenance.", "All analytical claims use the selected cohort."));
  }

  if (slides.length !== spec.slides) throw new Error(`${variant}: expected ${spec.slides} slides, built ${slides.length}`);
  return slides;
}

function visibleTextboxText(ndjson) {
  return String(ndjson || "")
    .split(/\r?\n/)
    .filter(Boolean)
    .flatMap((line) => {
      try {
        const record = JSON.parse(line);
        return record.kind === "textbox" && typeof record.text === "string" ? [record.text] : [];
      } catch {
        return [];
      }
    })
    .join("\n");
}

async function exportQa({ presentation, slides, workspace, candidatePptx, artifact, c }) {
  const previewDir = path.join(workspace, "preview", "final");
  const layoutDir = path.join(workspace, "layout", "final");
  await fs.mkdir(previewDir, { recursive: true });
  await fs.mkdir(layoutDir, { recursive: true });
  const previewPaths = [];
  for (let index = 0; index < slides.length; index += 1) {
    const padded = padSlideNumber(index + 1);
    const previewPath = path.join(previewDir, `slide-${padded}.png`);
    await saveBlobToFile(await presentation.export({ slide: slides[index], format: "png", scale: 2 }), previewPath);
    await saveBlobToFile(await presentation.export({ slide: slides[index], format: "layout" }), path.join(layoutDir, `slide-${padded}.layout.json`));
    previewPaths.push(previewPath);
  }
  await saveBlobToFile(await presentation.export({ format: "webp", montage: true, scale: 1 }), path.join(workspace, "final-montage.webp"));
  const inspection = await presentation.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 500000 });
  await writeTextAtomic(path.join(workspace, "final-inspect.ndjson"), inspection.ndjson || "");
  const visibleText = visibleTextboxText(inspection.ndjson);
  const visibleAnchors = [
    String(c.episodes),
    String(c.residents),
    `${c.uti} operational UTI`,
    `${c.notUti} operational Not_UTI`,
    String(c.directPairs),
    String(c.transitions),
    String(c.transitionsAtThreshold),
    `${c.toUtiAtThreshold} of ${c.toUti}`,
    `${c.nearMiss}-row near-miss`,
    "RQ01",
    "RQ10",
  ];
  const missingAnchors = visibleAnchors.filter((anchor) => !visibleText.includes(anchor));
  if (missingAnchors.length) {
    throw new Error(`Visible claim-anchor check failed: ${missingAnchors.join(", ")}`);
  }
  const pptx = await artifact.PresentationFile.exportPptx(presentation);
  await pptx.save(candidatePptx);
  await assertNoEmptyPlaceholders(candidatePptx);
  await assertPptxHasNoForbiddenContent(candidatePptx);
  await run(PYTHON, [path.join(SKILL_DIR, "container_tools", "slides_test.py"), candidatePptx], { echo: true });
  await run(PYTHON, [path.join(SKILL_DIR, "template_following_scripts", "make_contact_sheet.py"), "--output", path.join(workspace, "final-contact-sheet.png"), ...previewPaths]);
  await writeJsonAtomic(path.join(workspace, "qa-manifest.json"), {
    slide_count: slides.length,
    rendered_slide_previews: previewPaths,
    layout_exports: previewPaths.map((_, index) => path.join(layoutDir, `slide-${padSlideNumber(index + 1)}.layout.json`)),
    visible_claim_anchors: visibleAnchors,
    programmatic_checks: {
      empty_placeholders: "PASS",
      forbidden_content: "PASS",
      slides_test: "PASS",
      visible_claim_anchors: "PASS",
    },
    required_visual_review: "Inspect every full-size rendered slide and record/fix any unintended overlap, clipping, wrapping, or low-legibility chart text before release.",
  });
  return { previewDir, layoutDir };
}

export async function buildReviewDeck(options = {}) {
  const variant = requireVariant(options.variant);
  const spec = VARIANTS[variant];
  const projectRoot = path.resolve(options.projectRoot || DEFAULT_PROJECT_ROOT);
  const devFixture = Boolean(options.devFixture);
  const workspace = path.resolve(options.workspace || path.join(os.tmpdir(), "ruti-longcycler-review-decks", variant));
  const canonicalOutput = path.join(CANONICAL_ROOT, spec.family, "output", spec.output);
  const outputPptx = path.resolve(options.output || canonicalOutput);
  if (devFixture && outputPptx.startsWith(CANONICAL_ROOT + path.sep)) {
    throw new Error("Development fixtures may only write outside the canonical presentation tree.");
  }
  const { registry, registryPath, registrySha256 } = await loadRegistry({ projectRoot, registryPath: options.registry, devFixture });
  const c = releaseCounts(registry);
  const plots = await resolvePlots(registry, projectRoot);

  await fs.rm(workspace, { recursive: true, force: true });
  await ensureArtifactToolWorkspace(workspace);
  const artifact = await importArtifactTool(workspace);
  const presentation = artifact.Presentation.create({ slideSize: SIZE });
  const slides = await buildNarrative({ presentation, artifact, workspace, variant, c, registry, registryPath, registrySha: registrySha256, plots });
  const candidatePptx = path.join(workspace, "candidate.pptx");
  const qa = await exportQa({ presentation, slides, workspace, candidatePptx, artifact, c });
  await copyFileAtomic(candidatePptx, outputPptx);
  await writeJsonAtomic(path.join(workspace, "generation-report.json"), {
    variant,
    output: outputPptx,
    slide_count: slides.length,
    registry: registryPath,
    registry_sha256: registrySha256,
    analytical_cohort: registry.analytical_cohort,
    direct_pairs: registry.direct_pairs,
    adjacent_transitions: registry.adjacent_transitions,
    mechanism_casebook: registry.mechanism_casebook,
    near_miss_audit: registry.near_miss_audit,
    research_questions: registry.research_questions,
    method_contract: registry.method_contract,
    source_plot_files: plots,
  });
  return { variant, outputPptx, workspace, candidatePptx, previewDir: qa.previewDir, layoutDir: qa.layoutDir, slideCount: slides.length };
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const result = await buildReviewDeck({
    variant: args.variant,
    projectRoot: args["project-root"],
    registry: args.registry,
    workspace: args.workspace,
    output: args.output,
    devFixture: args["dev-fixture"],
  });
  console.log(JSON.stringify(result, null, 2));
}

if (import.meta.url === `file://${process.argv[1]}`) {
  main().catch((error) => {
    console.error(error.stack || error.message || String(error));
    process.exit(1);
  });
}

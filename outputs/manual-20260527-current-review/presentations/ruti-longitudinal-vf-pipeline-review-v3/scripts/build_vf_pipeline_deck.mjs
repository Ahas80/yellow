#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { spawnSync } from "node:child_process";

import {
  createSlideContext,
  ensureArtifactToolWorkspace,
  importArtifactTool,
  padSlideNumber,
  saveBlobToFile,
} from "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.521.10419/skills/presentations/scripts/artifact_tool_utils.mjs";

const WORKSPACE = "/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v3";
const STARTER = path.join(WORKSPACE, "template-starter.pptx");
const OUTPUT = path.join(WORKSPACE, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx");
const PREVIEW_DIR = path.join(WORKSPACE, "preview", "final");
const LAYOUT_DIR = path.join(WORKSPACE, "layout", "final");
const CONTACT_SHEET = path.join(WORKSPACE, "preview", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_contact_sheet.png");
const SKILL_DIR = "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.521.10419/skills/presentations";

const ASSETS = {
  vfTopGenes: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_top_gene_prevalence.png",
  moduleGeneCounts: "/Users/Aamir/Desktop/rUTIs/plots/vf/module_gene_counts.png",
  withinHostJaccard: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_within_host_jaccard_distribution.png",
  geneGainLoss: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_gene_gain_loss_consecutive_pairs.png",
  jaccardSameSt: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_jaccard_same_vs_different_st.png",
  geneModelEvidence: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_gene_screening_vs_model_evidence.png",
  modulePrevalence: "/Users/Aamir/Desktop/rUTIs/plots/vf/module_prevalence_by_status.png",
  main1: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Main_Figure_1_primary_denominator_and_uncertainty.png",
  main2: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Main_Figure_2_not_uti_to_uti_mechanism_casebook.png",
  main3: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Main_Figure_3_strain_stability_and_host_context.png",
  main4: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Main_Figure_4_global_vf_signal_and_robustness.png",
  pcoaSt: "/Users/Aamir/Desktop/rUTIs/plots/vf/vf_pcoa_jaccard_ST.png",
  suppS2: "/Users/Aamir/Desktop/rUTIs/plots/final_figures/Supplementary_Figure_S2_accessory_plasmid_amr_changes.png",
};

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
  red: "#B91C1C",
};

function slidesFromPresentation(presentation) {
  if (Array.isArray(presentation.slides?.items)) return presentation.slides.items;
  return Array.from({ length: presentation.slides.count }, (_, index) => presentation.slides.getItem(index));
}

function deleteCollection(collection) {
  const items = Array.isArray(collection?.items) ? [...collection.items] : [];
  for (const item of items) {
    if (typeof item.delete === "function") item.delete();
    else if (item.id && typeof collection.deleteById === "function") collection.deleteById(item.id);
  }
}

function clearSlide(slide) {
  if (slide.shapes?.deleteAll) slide.shapes.deleteAll();
  deleteCollection(slide.images);
  deleteCollection(slide.tables);
  deleteCollection(slide.charts);
  if (slide.speakerNotes?.clear) slide.speakerNotes.clear();
}

function base(ctx, slide, eyebrow, title, subtitle, opts = {}) {
  ctx.addShape(slide, { x: 0, y: 0, w: ctx.W, h: ctx.H, fill: "#FFFFFF", line: ctx.line() });
  ctx.addShape(slide, { x: 0, y: 0, w: 12, h: ctx.H, fill: opts.rail || C.blue, line: ctx.line() });
  ctx.addText(slide, { x: 54, y: 34, w: 900, h: 20, text: eyebrow.toUpperCase(), fontSize: 12, bold: true, color: C.muted, typeface: "Aptos" });
  ctx.addText(slide, { x: 54, y: 58, w: 1120, h: 44, text: title, fontSize: opts.titleSize || 30, bold: true, color: C.ink, typeface: "Aptos Display" });
  if (subtitle) ctx.addText(slide, { x: 54, y: 104, w: 1120, h: 28, text: subtitle, fontSize: 16, color: C.muted, typeface: "Aptos" });
  ctx.addShape(slide, { x: 54, y: 140, w: 1172, h: 1.4, fill: C.line, line: ctx.line() });
  ctx.addShape(slide, { x: 54, y: 674, w: 1172, h: 1.2, fill: C.line, line: ctx.line() });
  const footer = opts.footer || "Current pipeline outputs | VF analysis of longitudinal urinary E. coli isolates";
  ctx.addText(slide, { x: 54, y: 682, w: 1000, h: 20, text: footer, fontSize: 10.5, color: C.muted });
  ctx.addText(slide, { x: 1160, y: 682, w: 66, h: 20, text: String(ctx.slideNumber).padStart(2, "0"), fontSize: 10.5, color: C.muted, align: "right" });
}

function label(ctx, slide, x, y, text, fill = C.slateFill, color = C.ink, w = 148) {
  ctx.addShape(slide, { x, y, w, h: 26, fill, line: { style: "solid", fill: fill, width: 0 } });
  ctx.addText(slide, { x: x + 10, y: y + 5, w: w - 20, h: 16, text, fontSize: 10.5, bold: true, color, align: "center" });
}

function card(ctx, slide, x, y, w, h, title, body, opts = {}) {
  ctx.addShape(slide, {
    x, y, w, h,
    fill: opts.fill || C.light,
    line: { style: "solid", fill: opts.line || C.line, width: 1 },
    geometry: "roundRect",
  });
  if (opts.kicker) ctx.addText(slide, { x: x + 18, y: y + 14, w: w - 36, h: 16, text: opts.kicker, fontSize: 10.5, bold: true, color: opts.kickerColor || C.muted });
  const titleY = y + (opts.kicker ? 34 : 18);
  const titleH = opts.titleH || 24;
  const bodyY = opts.bodyY || y + (opts.kicker ? 62 : 48);
  ctx.addText(slide, { x: x + 18, y: titleY, w: w - 36, h: titleH, text: title, fontSize: opts.titleSize || 18, bold: true, color: opts.titleColor || C.ink });
  ctx.addText(slide, { x: x + 18, y: bodyY, w: w - 36, h: Math.max(10, y + h - bodyY - 14), text: body, fontSize: opts.bodySize || 13.5, color: opts.bodyColor || C.muted, insets: { left: 0, right: 0, top: 0, bottom: 0 } });
}

function metric(ctx, slide, x, y, w, value, labelText, fill, color = C.ink) {
  ctx.addShape(slide, { x, y, w, h: 86, fill, line: { style: "solid", fill, width: 0 }, geometry: "roundRect" });
  ctx.addText(slide, { x: x + 16, y: y + 14, w: w - 32, h: 32, text: value, fontSize: 28, bold: true, color, align: "center" });
  ctx.addText(slide, { x: x + 16, y: y + 52, w: w - 32, h: 20, text: labelText, fontSize: 11.5, color, align: "center" });
}

function arrow(ctx, slide, x1, y1, x2, y2, color = C.line) {
  const dx = x2 - x1;
  const dy = y2 - y1;
  const len = Math.max(1, Math.sqrt(dx * dx + dy * dy));
  const ang = Math.atan2(dy, dx);
  const cx = (x1 + x2) / 2;
  const cy = (y1 + y2) / 2;
  const line = ctx.addShape(slide, { x: cx - len / 2, y: cy - 1, w: len, h: 2, fill: color, line: ctx.line() });
  line.position.rotation = (ang * 180) / Math.PI;
  const head = ctx.addShape(slide, { x: x2 - 7, y: y2 - 5, w: 10, h: 10, fill: color, line: ctx.line(), geometry: "triangle" });
  head.position.rotation = (ang * 180) / Math.PI + 90;
}

async function addFigure(ctx, slide, imagePath, opts = {}) {
  const x = opts.x ?? 64;
  const y = opts.y ?? 158;
  const w = opts.w ?? 1040;
  const h = opts.h ?? 480;
  ctx.addShape(slide, { x: x - 8, y: y - 8, w: w + 16, h: h + 16, fill: "#FFFFFF", line: { style: "solid", fill: "#E5E7EB", width: 1 } });
  await ctx.addImage(slide, { x, y, w, h, path: imagePath, fit: "contain", alt: opts.alt || path.basename(imagePath) });
}

function slide1(ctx, slide) {
  ctx.addShape(slide, { x: 0, y: 0, w: ctx.W, h: ctx.H, fill: "#FFFFFF", line: ctx.line() });
  ctx.addShape(slide, { x: 0, y: 0, w: 16, h: ctx.H, fill: C.blue, line: ctx.line() });
  ctx.addText(slide, { x: 62, y: 58, w: 1050, h: 24, text: "CURRENT PIPELINE REVIEW  /  LONGITUDINAL VF ANALYSIS", fontSize: 12.5, bold: true, color: C.muted });
  ctx.addText(slide, { x: 62, y: 96, w: 1050, h: 92, text: "Virulence-factor profiling of longitudinal urinary E. coli isolates", fontSize: 39, bold: true, color: C.ink, typeface: "Aptos Display" });
  ctx.addText(slide, { x: 64, y: 190, w: 900, h: 30, text: "Current pipeline, descriptive longitudinal findings, and cautious clinical-status overlay", fontSize: 18, color: C.muted });
  card(ctx, slide, 64, 286, 340, 142, "What VF repertoire is present?", "Describe the observed virulence-factor gene space and curated modules in 556 VF/WGS-ready urinary E. coli isolates.", { fill: C.blueFill, line: "#B9D5F0", kicker: "Question 1", kickerColor: C.blue });
  card(ctx, slide, 462, 286, 340, 142, "How stable are profiles over time?", "Use consecutive within-resident isolate-pair comparisons to quantify VF similarity and gene gain/loss.", { fill: C.greenFill, line: "#BCE4C9", kicker: "Question 2", kickerColor: C.green });
  card(ctx, slide, 860, 286, 340, 142, "Clinical status overlay", "Treat UTI versus Not_UTI as a sparse, exploratory phenotype annotation rather than the whole study frame.", { fill: C.orangeFill, line: "#F4C98A", kicker: "Question 3", kickerColor: C.orange });
  metric(ctx, slide, 190, 500, 185, "556", "VF-ready episodes", C.blueFill, C.blue);
  metric(ctx, slide, 415, 500, 185, "227", "VF gene columns", C.slateFill, C.ink);
  metric(ctx, slide, 640, 500, 185, "394", "longitudinal pairs", C.greenFill, C.green);
  metric(ctx, slide, 865, 500, 185, "17", "UTI annotations", C.orangeFill, C.orange);
  ctx.addText(slide, { x: 62, y: 682, w: 900, h: 20, text: "Designed for a 22 min scientific review plus discussion", fontSize: 11, color: C.muted });
}

function slide2(ctx, slide) {
  base(ctx, slide, "Why this framing", "Longitudinal urinary E. coli VF analysis asks about repertoire and stability", "Clinical status is an annotation layer; the core object is repeated isolate VF profiles.");
  const y = 220;
  card(ctx, slide, 78, 184, 292, 310, "Repeated resident isolates", "Urinary E. coli can recur or persist across serial sampling. A one-time status comparison cannot separate stable carriage, replacement, and longitudinal change.", { fill: C.blueFill, line: "#BDD7EE", bodySize: 15 });
  ctx.addText(slide, { x: 430, y: 186, w: 400, h: 24, text: "Three observable longitudinal patterns", fontSize: 19, bold: true, color: C.ink });
  const items = [
    ["Persistent lineage", "High VF similarity across repeated isolates", C.greenFill, C.green],
    ["Replacement", "Different lineage or large genomic distance", C.rustFill, C.rust],
    ["Gene-content change", "Detected VF gains/losses between visits", C.slateFill, C.ink],
  ];
  items.forEach((d, i) => {
    const yy = y + i * 86;
    ctx.addShape(slide, { x: 430, y: yy, w: 290, h: 54, fill: d[2], line: { style: "solid", fill: d[2], width: 0 }, geometry: "roundRect" });
    ctx.addText(slide, { x: 448, y: yy + 10, w: 250, h: 18, text: d[0], fontSize: 15, bold: true, color: d[3] });
    ctx.addText(slide, { x: 448, y: yy + 31, w: 250, h: 16, text: d[1], fontSize: 11.5, color: C.muted });
  });
  card(ctx, slide, 790, 184, 374, 310, "Clinical overlay comes later", "UTI_Status still matters, but the sparse UTI denominator means clinical contrasts should be labelled exploratory. The VF pipeline first establishes what was measured, how profiles are curated, and how stable those profiles are within residents.", { fill: "#FFFFFF", line: C.line, titleColor: C.orange, bodySize: 15 });
  label(ctx, slide, 790, 526, "Not_UTI", C.blueFill, C.blue);
  label(ctx, slide, 970, 526, "UTI", C.orangeFill, C.orange);
  arrow(ctx, slide, 936, 539, 966, 539, C.line);
}

function slide3(ctx, slide) {
  base(ctx, slide, "Data asset", "VF-ready data link repeated urinary E. coli genomes to clinical episodes", "Dataset-level counts are VF-ready denominators; no disputed clinical participant total is displayed.");
  metric(ctx, slide, 80, 180, 190, "556", "episode-level selected assemblies", C.blueFill, C.blue);
  metric(ctx, slide, 300, 180, 190, "162", "participants in VF-ready dataset", C.slateFill, C.ink);
  metric(ctx, slide, 520, 180, 190, "227", "binary VF gene columns", C.greenFill, C.green);
  metric(ctx, slide, 740, 180, 190, "32", "curated modules", C.orangeFill, C.orange);
  metric(ctx, slide, 960, 180, 190, "18", "UPEC-candidate modules", C.rustFill, C.rust);
  ctx.addText(slide, { x: 98, y: 340, w: 300, h: 22, text: "Repeated isolate design", fontSize: 20, bold: true, color: C.ink });
  const xs = [150, 370, 590, 810, 1030];
  xs.forEach((x, i) => {
    ctx.addShape(slide, { x, y: 440, w: 70, h: 70, fill: i === 3 ? C.orangeFill : C.blueFill, line: { style: "solid", fill: i === 3 ? "#F4B15F" : "#9AC1E6", width: 2 }, geometry: "ellipse" });
    ctx.addText(slide, { x: x - 8, y: 468, w: 86, h: 16, text: i === 3 ? "UTI" : "Not_UTI", fontSize: 11, bold: true, color: i === 3 ? C.orange : C.blue, align: "center" });
    if (i < xs.length - 1) arrow(ctx, slide, x + 80, 475, xs[i + 1] - 10, 475, "#94A3B8");
    ctx.addText(slide, { x: x - 16, y: 525, w: 102, h: 16, text: `Episode ${i + 1}`, fontSize: 11.5, color: C.muted, align: "center" });
  });
  card(ctx, slide, 98, 558, 1050, 82, "Interpretation stance", "The unit of VF measurement is an episode-level isolate profile. Status labels can be overlaid, but repeated measures, lineage structure, and sparse UTI counts limit causal status claims.", { fill: C.light, bodySize: 12.5, titleSize: 14 });
}

function slide4(ctx, slide) {
  base(ctx, slide, "Current pipeline", "Numbered clinical-to-VF pipeline with visible denominators at each handoff", "Counts come from current denominator summaries and validation checks; assembler alternatives are labelled explicitly.", { titleSize: 28 });
  const stages = [
    ["1", "Clinical classification", "00a / 00b", "585 classified episodes before manual exclusions"],
    ["2", "Primary inclusion", "status_map.csv", "583 included: 18 UTI, 565 Not_UTI; 2 excluded"],
    ["3", "Assembly QC", "12a_wgs_qc.R", "1,291 assembly-level QC records incl. assembler alternatives"],
    ["4", "Canonical selection", "canonical_assembly_selection.csv", "556 selected episode-level assemblies"],
    ["5", "VF P/A matrix", "02_gene_presence_analysis.R", "556 episode-level VF rows"],
    ["6", "VF/model-ready", "22_vf_build_analysis_dataset.R", "556 episodes: 17 UTI, 539 Not_UTI; 27 clinical episodes lack VF-ready evidence"],
    ["7", "Feature framework", "26_vf_define_gene_modules.R", "227 VF gene columns -> 32 modules, incl. 18 UPEC-candidate modules"],
    ["8", "Longitudinal pairs", "24_vf_longitudinal_dynamics.R", "394 consecutive comparisons from 144 participants"],
    ["9", "Clinical application", "casebook outputs", "11 Not_UTI -> UTI transitions; 10 WGS/VF-linked; 1 missing endpoint"],
  ];
  stages.forEach((s, i) => {
    const col = i % 3;
    const row = Math.floor(i / 3);
    const x = 70 + col * 392;
    const y = 166 + row * 158;
    ctx.addShape(slide, { x, y, w: 344, h: 112, fill: i < 2 ? C.orangeFill : i < 6 ? C.blueFill : i === 7 ? C.greenFill : C.slateFill, line: { style: "solid", fill: "#D8DEE8", width: 1 }, geometry: "roundRect" });
    ctx.addShape(slide, { x: x + 16, y: y + 18, w: 30, h: 30, fill: C.ink, line: ctx.line(), geometry: "ellipse" });
    ctx.addText(slide, { x: x + 16, y: y + 24, w: 30, h: 16, text: s[0], fontSize: 11.5, bold: true, color: "#FFFFFF", align: "center" });
    ctx.addText(slide, { x: x + 58, y: y + 15, w: 260, h: 20, text: s[1], fontSize: 15, bold: true, color: C.ink });
    ctx.addText(slide, { x: x + 58, y: y + 38, w: 260, h: 15, text: s[2], fontSize: 10.5, color: C.muted });
    ctx.addText(slide, { x: x + 18, y: y + 64, w: 308, h: 35, text: s[3], fontSize: 11.5, color: C.ink });
    if (col < 2) arrow(ctx, slide, x + 348, y + 56, x + 382, y + 56, "#94A3B8");
  });
}

function slide5(ctx, slide) {
  base(ctx, slide, "VF representation", "VF features: genes, curated modules, and longitudinal similarity", "Module and score outputs help interpretation, but they are not validated disease-causality scores.");
  card(ctx, slide, 74, 184, 280, 140, "Binary VF gene matrix", "227 presence/absence gene columns per episode-level isolate.", { fill: C.blueFill, line: "#B9D5F0", titleColor: C.blue });
  card(ctx, slide, 500, 184, 280, 140, "Curated modules", "32 biological curation units, including 18 UPEC-candidate modules.", { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green });
  card(ctx, slide, 926, 184, 280, 140, "Outputs", "Burden, prevalence, Jaccard similarity, gain/loss, and exploratory models.", { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange });
  arrow(ctx, slide, 368, 252, 486, 252, "#94A3B8");
  arrow(ctx, slide, 794, 252, 912, 252, "#94A3B8");
  const matrixX = 110, matrixY = 398;
  for (let r = 0; r < 6; r++) {
    for (let c = 0; c < 10; c++) {
      const on = (r + c * 2) % 3 !== 0;
      ctx.addShape(slide, { x: matrixX + c * 19, y: matrixY + r * 17, w: 14, h: 12, fill: on ? C.blue : "#E2E8F0", line: ctx.line() });
    }
  }
  ctx.addText(slide, { x: 74, y: 518, w: 290, h: 36, text: "Rows are episode-level isolates; columns are VF genes detected from current selected assemblies.", fontSize: 12.5, color: C.muted });
  ["Adhesion", "Iron acquisition", "Secretion", "Toxin", "Capsule", "Unassigned"].forEach((name, i) => {
    ctx.addShape(slide, { x: 508 + (i % 2) * 142, y: 386 + Math.floor(i / 2) * 44, w: 122, h: 28, fill: i === 5 ? C.slateFill : C.greenFill, line: { style: "solid", fill: "#D8DEE8", width: 1 }, geometry: "roundRect" });
    ctx.addText(slide, { x: 516 + (i % 2) * 142, y: 393 + Math.floor(i / 2) * 44, w: 106, h: 14, text: name, fontSize: 10.5, color: i === 5 ? C.muted : C.green, align: "center", bold: true });
  });
  card(ctx, slide, 902, 380, 300, 146, "Interpretation caution", "Unassigned genes are 25.1% of the VF matrix. Interpret total burden separately from curated and UPEC-candidate summaries.", { fill: "#FFF7ED", line: "#FED7AA", titleColor: C.orange, bodySize: 13.5 });
  label(ctx, slide, 914, 548, "Descriptive framework", C.slateFill, C.muted, 200);
}

async function slideFigure(ctx, slide, slideSpec) {
  base(ctx, slide, slideSpec.eyebrow, slideSpec.title, slideSpec.subtitle, { titleSize: slideSpec.titleSize || 28, footer: slideSpec.footer });
  await addFigure(ctx, slide, slideSpec.image, { x: slideSpec.x || 62, y: slideSpec.y || 158, w: slideSpec.w || 930, h: slideSpec.h || 480, alt: slideSpec.alt });
  if (slideSpec.callout) card(ctx, slide, slideSpec.callout.x, slideSpec.callout.y, slideSpec.callout.w, slideSpec.callout.h, slideSpec.callout.title, slideSpec.callout.body, slideSpec.callout);
  if (slideSpec.evidence) label(ctx, slide, slideSpec.evidence.x || 1010, slideSpec.evidence.y || 160, slideSpec.evidence.text, slideSpec.evidence.fill || C.slateFill, slideSpec.evidence.color || C.muted, slideSpec.evidence.w || 178);
}

function slide12(ctx, slide) {
  base(ctx, slide, "Clinical annotation application 2", "Participant 20026 illustrates stable VF profile despite symptom emergence", "This is an applied phenotype example, not proof of host-state mechanism.");
  ctx.addShape(slide, { x: 120, y: 312, w: 960, h: 4, fill: "#94A3B8", line: ctx.line() });
  ctx.addShape(slide, { x: 196, y: 268, w: 154, h: 96, fill: C.blueFill, line: { style: "solid", fill: "#9AC1E6", width: 2 }, geometry: "roundRect" });
  ctx.addText(slide, { x: 214, y: 286, w: 118, h: 20, text: "T3", fontSize: 17, bold: true, color: C.blue, align: "center" });
  ctx.addText(slide, { x: 214, y: 312, w: 118, h: 20, text: "Not_UTI", fontSize: 14, bold: true, color: C.blue, align: "center" });
  ctx.addText(slide, { x: 206, y: 338, w: 132, h: 16, text: "culture-supported bacteriuria", fontSize: 9.5, color: C.muted, align: "center" });
  ctx.addShape(slide, { x: 858, y: 268, w: 154, h: 96, fill: C.orangeFill, line: { style: "solid", fill: "#F4B15F", width: 2 }, geometry: "roundRect" });
  ctx.addText(slide, { x: 876, y: 286, w: 118, h: 20, text: "UTI-1", fontSize: 17, bold: true, color: C.orange, align: "center" });
  ctx.addText(slide, { x: 876, y: 312, w: 118, h: 20, text: "UTI", fontSize: 14, bold: true, color: C.orange, align: "center" });
  ctx.addText(slide, { x: 868, y: 338, w: 132, h: 16, text: "symptom emergence", fontSize: 9.5, color: C.muted, align: "center" });
  arrow(ctx, slide, 360, 314, 848, 314, "#94A3B8");
  label(ctx, slide, 505, 278, "42 days", C.slateFill, C.ink, 120);
  metric(ctx, slide, 210, 432, 160, "5", "SNPs", C.greenFill, C.green);
  metric(ctx, slide, 410, 432, 160, "+0/-0", "VF genes", C.greenFill, C.green);
  metric(ctx, slide, 610, 432, 160, "stable", "VF modules", C.greenFill, C.green);
  metric(ctx, slide, 810, 432, 160, "strong", "same-strain evidence", C.greenFill, C.green);
  card(ctx, slide, 80, 552, 520, 78, "Transition casebook context", "11 clinical Not_UTI -> UTI transitions; 10 WGS/VF-linked; 1 missing endpoint.", { fill: C.light, titleSize: 14, bodySize: 12 });
  card(ctx, slide, 650, 552, 520, 78, "Mechanism summary", "4 stable-profile transitions and 3 replacement-consistent transitions; buckets organise evidence, not mechanism proof.", { fill: C.light, titleSize: 14, bodySize: 12 });
}

function slide13(ctx, slide) {
  base(ctx, slide, "Close", "What this VF-first review establishes and what remains uncertain", "Use these as discussion prompts and handover priorities.");
  card(ctx, slide, 82, 184, 330, 296, "Evidence", "The current pipeline produces a coherent VF-ready longitudinal dataset: 556 isolate profiles, 227 VF gene columns, 32 modules, and 394 consecutive within-resident comparisons.", { fill: C.blueFill, line: "#B9D5F0", titleColor: C.blue, bodySize: 15 });
  card(ctx, slide, 474, 184, 330, 296, "Boundary", "VF profiles are often stable, but presence/absence does not measure expression or activity. Clinical status associations are exploratory because UTI counts are sparse.", { fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 15 });
  card(ctx, slide, 866, 184, 330, 296, "Next steps", "Use the pipeline outputs to prioritise lineage-aware longitudinal follow-up, expression/regulation hypotheses, and targeted review of transition case studies.", { fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 15 });
  ctx.addText(slide, { x: 100, y: 540, w: 1080, h: 26, text: "Main discussion question: what additional evidence would distinguish stable carriage plus host-state change from unmeasured bacterial regulation?", fontSize: 18, bold: true, color: C.ink, align: "center" });
}

function slide14(ctx, slide) {
  base(ctx, slide, "Appendix / handover", "Practical map: where to enter and rerun the VF pipeline", "Keep this as the operational handover slide for the colleague.");
  const rows = [
    ["Clinical key", "results/clinical/status_map.csv", "Primary UTI_Status annotation and audit fields"],
    ["VF matrix", "results/vf/vf_pa_all.csv", "Episode-level gene presence/absence rows"],
    ["Analysis dataset", "results/vf/vf_analysis_ready.csv", "Canonical 556-row VF/model-ready table"],
    ["Modules", "results/vf/gene_module_map.csv", "Gene-to-module curation and unassigned review"],
    ["Longitudinal", "results/vf/vf_longitudinal_transitions.csv", "Consecutive within-resident VF comparisons"],
    ["Interpretation", "results/vf/vf_longitudinal_summary.txt", "Stable-profile and gain/loss summaries"],
    ["Final figures", "35_final_figure_pack.R", "Entry point for current final figure pack"],
  ];
  rows.forEach((r, i) => {
    const y = 164 + i * 64;
    ctx.addShape(slide, { x: 76, y, w: 1070, h: 48, fill: i % 2 ? "#FFFFFF" : C.light, line: { style: "solid", fill: "#E5E7EB", width: 1 } });
    ctx.addText(slide, { x: 96, y: y + 9, w: 150, h: 18, text: r[0], fontSize: 13, bold: true, color: C.blue });
    ctx.addText(slide, { x: 270, y: y + 9, w: 350, h: 18, text: r[1], fontSize: 11.5, color: C.ink, typeface: "Aptos Mono" });
    ctx.addText(slide, { x: 650, y: y + 9, w: 470, h: 28, text: r[2], fontSize: 12.5, color: C.muted });
  });
}

function slide22(ctx, slide) {
  base(ctx, slide, "Appendix / evidence registry", "Sources of truth and lookup resources for questions", "Open these files when the discussion needs provenance rather than a new analysis.");
  const left = [
    ["Counts and validation", "results/final_figures/final_figure_validation_checks.csv"],
    ["Denominator flow", "results/qc/pipeline_denominator_summary.csv"],
    ["VF diagnostics", "results/vf/vf_dataset_diagnostics.txt"],
    ["VF figure registry", "results/vf/vf_figure_index.csv"],
    ["VF visual audit", "results/vf/vf_visualisation_audit.csv"],
  ];
  const right = [
    ["Longitudinal summary", "results/vf/vf_longitudinal_summary.txt"],
    ["Module notes", "results/vf/vf_module_definition_notes.md"],
    ["Casebook", "results/mechanism/not_uti_to_uti_casebook.md"],
    ["Pipeline explanation", "docs/pipeline_architecture.md"],
    ["Dense backup lookup", "results/longitudinal/swimmer_plot.png"],
  ];
  const renderList = (items, x) => items.forEach((r, i) => {
    const y = 174 + i * 72;
    ctx.addShape(slide, { x, y, w: 500, h: 54, fill: C.light, line: { style: "solid", fill: "#E5E7EB", width: 1 }, geometry: "roundRect" });
    ctx.addText(slide, { x: x + 18, y: y + 10, w: 190, h: 16, text: r[0], fontSize: 12.5, bold: true, color: C.ink });
    ctx.addText(slide, { x: x + 18, y: y + 30, w: 455, h: 14, text: r[1], fontSize: 10.5, color: C.muted, typeface: "Aptos Mono" });
  });
  renderList(left, 82);
  renderList(right, 680);
  card(ctx, slide, 190, 570, 900, 56, "Guardrail for Q&A", "Use Supplementary Figures S1 and S3-S5 for transition/sensitivity questions; do not convert backup diagnostics into confirmatory VF, AMR, or mechanism claims.", { fill: C.orangeFill, line: "#FED7AA", titleSize: 14, bodySize: 12 });
}

const figureSlides = {
  6: { eyebrow: "VF repertoire", title: "Common VF genes describe the observed repertoire, not disease association", subtitle: "Top-gene prevalence across 556 VF/WGS-linked urinary E. coli isolates.", image: ASSETS.vfTopGenes, w: 900, h: 472, callout: { x: 1010, y: 186, w: 180, h: 170, title: "Read as", body: "Descriptive prevalence ranking only. Gene presence does not imply expression, function, or UTI association.", fill: C.blueFill, line: "#B9D5F0", titleColor: C.blue, bodySize: 12.5 }, evidence: { text: "Descriptive", fill: C.blueFill, color: C.blue } },
  7: { eyebrow: "VF framework", title: "Curated modules organise genes into interpretable biological systems", subtitle: "Modules are curation units for navigation; they are not validated disease-causality scores.", image: ASSETS.moduleGeneCounts, w: 870, h: 480, callout: { x: 990, y: 188, w: 205, h: 195, title: "Why it matters", body: "The framework separates adhesion, iron acquisition, secretion, capsule/surface, toxin and unassigned material before longitudinal or clinical overlays.", fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 12.5 }, evidence: { text: "Descriptive framework", fill: C.greenFill, color: C.green, w: 202 } },
  8: { eyebrow: "Longitudinal stability", title: "Most repeated within-resident VF profiles are highly stable", subtitle: "394 consecutive comparisons from 144 participants; median Jaccard similarity 1.000; 62.4% no VF change.", image: ASSETS.withinHostJaccard, w: 905, h: 470, callout: { x: 1005, y: 188, w: 190, h: 190, title: "Central finding", body: "The dominant pattern is VF profile conservation across repeated urinary E. coli isolates from the same resident.", fill: C.greenFill, line: "#BCE4C9", titleColor: C.green, bodySize: 12.5 }, evidence: { text: "Descriptive longitudinal", fill: C.greenFill, color: C.green, w: 210 } },
  9: { eyebrow: "Longitudinal change", title: "When profiles change, gain/loss summaries flag candidates for follow-up", subtitle: "Observed gain/loss may reflect replacement, assembly/calling differences, or true gene-content change.", image: ASSETS.geneGainLoss, x: 58, y: 156, w: 1030, h: 482, callout: { x: 1102, y: 186, w: 118, h: 218, title: "Caution", body: "Do not overread gain/loss without lineage and genome-distance context.", fill: C.slateFill, line: "#CBD5E1", titleColor: C.ink, bodySize: 11.5 }, evidence: { x: 1096, y: 420, text: "Descriptive", fill: C.slateFill, color: C.muted, w: 124 } },
  10: { eyebrow: "Lineage context", title: "Sequence-type consistency helps interpret VF stability", subtitle: "Same-ST comparisons support persistent-lineage context, but ST agreement alone does not prove same strain.", image: ASSETS.jaccardSameSt, w: 900, h: 470, callout: { x: 1005, y: 186, w: 190, h: 206, title: "Use as diagnostic", body: "Lineage context is essential before interpreting profile changes or clinical-status contrasts.", fill: C.slateFill, line: "#CBD5E1", titleColor: C.ink, bodySize: 12.5 }, evidence: { text: "Diagnostic / descriptive", fill: C.slateFill, color: C.muted, w: 220 } },
  11: { eyebrow: "Clinical annotation application 1", title: "UTI status overlay remains exploratory under the sparse denominator", subtitle: "556 model-ready episodes: 17 UTI and 539 Not_UTI; no global VF association remains significant after FDR correction.", image: ASSETS.geneModelEvidence, w: 890, h: 470, callout: { x: 1005, y: 186, w: 194, h: 210, title: "Safe claim", body: "Nominal screens and model diagnostics are hypothesis-generating. They do not establish a robust VF status association.", fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 12.5 }, evidence: { text: "Exploratory phenotype overlay", fill: C.orangeFill, color: C.orange, w: 230 } },
  15: { eyebrow: "Appendix", title: "Clinical phenotype denominator and definition", subtitle: "Clinical status is the annotation layer used for later exploratory overlays.", image: ASSETS.main1, x: 68, y: 158, w: 1080, h: 482, evidence: { x: 972, y: 160, text: "Denominator audit", fill: C.slateFill, color: C.muted, w: 200 } },
  16: { eyebrow: "Appendix", title: "Module prevalence by status is an exploratory annotation view", subtitle: "Read as descriptive status context; repeated measures and lineage are not fully resolved here.", image: ASSETS.modulePrevalence, x: 66, y: 158, w: 1040, h: 480, callout: { x: 1120, y: 190, w: 100, h: 165, title: "Label", body: "Exploratory clinical annotation.", fill: C.orangeFill, line: "#F4C98A", titleColor: C.orange, bodySize: 11 } },
  17: { eyebrow: "Appendix", title: "Full transition mechanism casebook", subtitle: "Available for questions about Not_UTI -> UTI transitions; buckets organise evidence, not proof.", image: ASSETS.main2, x: 62, y: 158, w: 1110, h: 482, evidence: { x: 980, y: 160, text: "Descriptive casebook", fill: C.slateFill, color: C.muted, w: 210 } },
  18: { eyebrow: "Appendix", title: "Stable strain and changing clinical state", subtitle: "Low SNP distance and stable VF profiles support, but do not prove, host-state or regulation hypotheses.", image: ASSETS.main3, x: 62, y: 158, w: 1110, h: 482, evidence: { x: 956, y: 160, text: "Transition context", fill: C.slateFill, color: C.muted, w: 210 } },
  19: { eyebrow: "Appendix", title: "Population-level robustness boundary", subtitle: "No adjusted VF association is confirmatory with only 17 VF-ready UTI episodes.", image: ASSETS.main4, x: 62, y: 158, w: 1110, h: 482, evidence: { x: 950, y: 160, text: "Exploratory diagnostics", fill: C.orangeFill, color: C.orange, w: 230 } },
  20: { eyebrow: "Appendix", title: "Lineage structure is an interpretation check", subtitle: "Jaccard PCoA of module profiles by sequence type flags lineage structure before status claims.", image: ASSETS.pcoaSt, x: 68, y: 158, w: 1020, h: 482, callout: { x: 1102, y: 190, w: 118, h: 190, title: "Use for", body: "Questions about ST and VF profile clustering.", fill: C.slateFill, line: "#CBD5E1", titleColor: C.ink, bodySize: 11.5 } },
  21: { eyebrow: "Appendix", title: "AMR backup only: exploratory transition-level ResFinder context", subtitle: "This is not a global VF-plus-AMR association analysis and should stay in Q&A backup.", image: ASSETS.suppS2, x: 62, y: 158, w: 1110, h: 482, evidence: { x: 890, y: 160, text: "AMR backup only", fill: C.rustFill, color: C.rust, w: 190 } },
};

async function build() {
  await ensureArtifactToolWorkspace(WORKSPACE);
  const { FileBlob, PresentationFile } = await importArtifactTool(WORKSPACE);
  const presentation = await PresentationFile.importPptx(await FileBlob.load(STARTER));
  const slides = slidesFromPresentation(presentation);
  if (slides.length !== 22) throw new Error(`Expected 22 slides, found ${slides.length}`);
  await fs.mkdir(PREVIEW_DIR, { recursive: true });
  await fs.mkdir(LAYOUT_DIR, { recursive: true });

  for (let i = 0; i < slides.length; i += 1) {
    const slide = slides[i];
    clearSlide(slide);
    const ctx = createSlideContext({ FileBlob, PresentationFile }, { slideSize: { width: 1280, height: 720 }, slideNumber: i + 1, workspaceDir: WORKSPACE, assetDir: path.join(WORKSPACE, "assets"), outputDir: path.dirname(OUTPUT) });
    if (i === 0) slide1(ctx, slide);
    else if (i === 1) slide2(ctx, slide);
    else if (i === 2) slide3(ctx, slide);
    else if (i === 3) slide4(ctx, slide);
    else if (i === 4) slide5(ctx, slide);
    else if (figureSlides[i + 1]) await slideFigure(ctx, slide, figureSlides[i + 1]);
    else if (i === 11) slide12(ctx, slide);
    else if (i === 12) slide13(ctx, slide);
    else if (i === 13) slide14(ctx, slide);
    else if (i === 21) slide22(ctx, slide);
    else throw new Error(`No renderer for slide ${i + 1}`);
    if (slide.speakerNotes?.setText) slide.speakerNotes.setText("");
  }

  const previewPaths = [];
  for (let i = 0; i < slides.length; i += 1) {
    const padded = padSlideNumber(i + 1);
    const preview = await presentation.export({ slide: slides[i], format: "png", scale: 1 });
    const previewPath = path.join(PREVIEW_DIR, `slide-${padded}.png`);
    await saveBlobToFile(preview, previewPath);
    previewPaths.push(previewPath);
    const layout = await presentation.export({ slide: slides[i], format: "layout" });
    await saveBlobToFile(layout, path.join(LAYOUT_DIR, `slide-${padded}.layout.json`));
  }

  const contact = spawnSync(
    process.env.PYTHON || "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3",
    [path.join(SKILL_DIR, "scripts", "make_contact_sheet.py"), "--output", CONTACT_SHEET, ...previewPaths],
    { encoding: "utf8" },
  );
  if (contact.status !== 0) throw new Error([contact.stdout, contact.stderr].filter(Boolean).join("\n"));

  await fs.mkdir(path.dirname(OUTPUT), { recursive: true });
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(OUTPUT);
  console.log(OUTPUT);
}

build().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

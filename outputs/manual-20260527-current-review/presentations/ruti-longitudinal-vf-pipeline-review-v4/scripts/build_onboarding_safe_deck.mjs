#!/usr/bin/env node

import { createRequire } from "module";

const require = createRequire(import.meta.url);
const fs = require("fs");
const fsp = require("fs/promises");
const path = require("path");
const { spawnSync } = require("child_process");
const sharp = require("/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp");
const pptxgen = require("/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/pptxgenjs");

const WORKSPACE = "/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v4";
const BASE_DIR = path.join(WORKSPACE, "preview", "base_shifted");
const FINAL_DIR = path.join(WORKSPACE, "preview", "final_25");
const OUTPUT = path.join(WORKSPACE, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx");
const CONTACT_SHEET = path.join(WORKSPACE, "preview", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_contact_sheet.png");
const SKILL_DIR = "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.521.10419/skills/presentations";
const PYTHON = "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";

const C = {
  ink: "#111827",
  muted: "#64748B",
  line: "#CBD5E1",
  light: "#F8FAFC",
  blue: "#2B6CB0",
  blueFill: "#E8F1FA",
  orange: "#D97706",
  orangeFill: "#FFF4E6",
  green: "#2F855A",
  greenFill: "#EAF7EF",
  rust: "#B65A3C",
  rustFill: "#FBEDE8",
  slateFill: "#EEF3F8",
};

function esc(text) {
  return String(text)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function wrapText(text, maxChars) {
  const words = String(text).split(/\s+/);
  const lines = [];
  let current = "";
  for (const word of words) {
    const candidate = current ? `${current} ${word}` : word;
    if (candidate.length > maxChars && current) {
      lines.push(current);
      current = word;
    } else {
      current = candidate;
    }
  }
  if (current) lines.push(current);
  return lines;
}

function textBlock(text, x, y, width, opts = {}) {
  const size = opts.size || 18;
  const lineHeight = opts.lineHeight || Math.round(size * 1.25);
  const color = opts.color || C.ink;
  const weight = opts.bold ? 700 : 400;
  const anchor = opts.anchor || "start";
  const lines = Array.isArray(text) ? text : wrapText(text, opts.maxChars || 38);
  const tspans = lines.map((line, i) => `<tspan x="${x}" dy="${i === 0 ? 0 : lineHeight}">${esc(line)}</tspan>`).join("");
  return `<text x="${x}" y="${y}" width="${width}" font-family="Aptos, Arial, sans-serif" font-size="${size}" font-weight="${weight}" fill="${color}" text-anchor="${anchor}">${tspans}</text>`;
}

function footer(slideNumber) {
  return `
    <line x1="54" y1="675" x2="1226" y2="675" stroke="${C.line}" stroke-width="1.2"/>
    ${textBlock("Current pipeline outputs | VF analysis of longitudinal urinary E. coli isolates", 54, 692, 900, { size: 10.5, color: C.muted, maxChars: 120 })}
    ${textBlock(String(slideNumber).padStart(2, "0"), 1226, 692, 40, { size: 10.5, color: C.muted, anchor: "end" })}
  `;
}

function base(eyebrow, title, subtitle, slideNumber) {
  return `
    <rect x="0" y="0" width="1280" height="720" fill="#FFFFFF"/>
    <rect x="0" y="0" width="12" height="720" fill="${C.blue}"/>
    ${textBlock(eyebrow.toUpperCase(), 54, 47, 900, { size: 12, color: C.muted, bold: true, maxChars: 120 })}
    ${textBlock(title, 54, 87, 1120, { size: 30, color: C.ink, bold: true, lineHeight: 36, maxChars: 66 })}
    ${textBlock(subtitle, 54, 120, 1120, { size: 16, color: C.muted, maxChars: 118 })}
    <line x1="54" y1="141" x2="1226" y2="141" stroke="${C.line}" stroke-width="1.4"/>
    ${footer(slideNumber)}
  `;
}

function roundedRect(x, y, w, h, fill, stroke = C.line, r = 10) {
  return `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${r}" fill="${fill}" stroke="${stroke}" stroke-width="1.2"/>`;
}

function arrow(x1, y1, x2, y2, color = "#94A3B8") {
  return `
    <line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="2.2"/>
    <polygon points="${x2},${y2} ${x2 - 9},${y2 - 5} ${x2 - 9},${y2 + 5}" fill="${color}"/>
  `;
}

function slide2Svg() {
  const nodes = [
    ["Clinical episode", "A urine-sampling moment with status, symptoms, culture context.", C.orangeFill, C.orange],
    ["Urinary E. coli isolate", "The bacterium selected from that episode for sequencing.", C.blueFill, C.blue],
    ["Genome assembly", "WGS reads are cleaned and assembled into a genome representation.", C.blueFill, C.blue],
    ["VF gene row", "Each isolate becomes detected / not detected across VF genes.", C.greenFill, C.green],
    ["Profile outputs", "Genes become modules, burden, similarity, and gain/loss summaries.", C.greenFill, C.green],
    ["Overlay & interpret", "Longitudinal patterns come first; clinical status is added cautiously.", C.orangeFill, C.orange],
  ];
  const boxW = 174;
  const gap = 24;
  const startX = 54;
  const y = 222;
  let body = base(
    "Project orientation",
    "The YELLOW routine in one plain-English loop",
    "One episode becomes one bacterial genome profile, then repeated profiles are compared over time.",
    2,
  );
  nodes.forEach((node, i) => {
    const x = startX + i * (boxW + gap);
    body += roundedRect(x, y, boxW, 168, node[2], node[3]);
    body += textBlock(String(i + 1), x + 16, y + 30, 30, { size: 14, bold: true, color: node[3] });
    body += textBlock(node[0], x + 18, y + 62, boxW - 36, { size: 17, bold: true, color: C.ink, lineHeight: 21, maxChars: 17 });
    body += textBlock(node[1], x + 18, y + 100, boxW - 36, { size: 11.5, color: C.muted, lineHeight: 15, maxChars: 23 });
    if (i < nodes.length - 1) body += arrow(x + boxW + 3, y + 84, x + boxW + gap - 5, y + 84);
  });
  body += roundedRect(92, 472, 318, 98, C.blueFill, "#B9D5F0");
  body += textBlock("What we measure", 112, 502, 280, { size: 16, bold: true, color: C.blue });
  body += textBlock("Presence/absence of virulence-factor genes in each selected urinary E. coli assembly.", 112, 530, 270, { size: 12.5, color: C.muted, maxChars: 42 });
  body += roundedRect(482, 472, 318, 98, C.greenFill, "#BCE4C9");
  body += textBlock("What we compare", 502, 502, 280, { size: 16, bold: true, color: C.green });
  body += textBlock("Repeated profiles within residents: how similar, what changed, and where replacement is plausible.", 502, 530, 270, { size: 12.5, color: C.muted, maxChars: 42 });
  body += roundedRect(872, 472, 318, 98, C.orangeFill, "#F4C98A");
  body += textBlock("Where status enters", 892, 502, 280, { size: 16, bold: true, color: C.orange });
  body += textBlock("UTI versus Not_UTI is an exploratory clinical annotation, not the whole analysis frame.", 892, 530, 270, { size: 12.5, color: C.muted, maxChars: 42 });
  return `<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="720">${body}</svg>`;
}

function slide3Svg() {
  const rows = [
    ["Clinical episodes", "583", "visits or events with clinical status and culture context", C.orangeFill, C.orange],
    ["VF-ready episodes", "556", "selected episode-level assemblies with VF-ready evidence", C.blueFill, C.blue],
    ["VF feature space", "227 -> 32", "VF gene columns grouped into curated modules", C.greenFill, C.green],
    ["Longitudinal pairs", "394", "consecutive within-resident comparisons", C.greenFill, C.green],
    ["Clinical transitions", "11", "Not_UTI -> UTI application cases", C.slateFill, C.ink],
  ];
  let body = base(
    "Project orientation",
    "What changes as data move through the pipeline",
    "Counts change because the unit changes: episode, assembly, gene/module, pair, or transition.",
    3,
  );
  body += roundedRect(72, 170, 1136, 74, "#FFFFFF", C.line);
  body += textBlock("Plain rule for reading the whole deck", 96, 202, 360, { size: 18, bold: true, color: C.ink });
  body += textBlock("Always ask: what is the row here? The denominator changes when we move from clinical records to genomes, VF features, paired comparisons, and transition cases.", 430, 198, 720, { size: 15, color: C.muted, maxChars: 92, lineHeight: 19 });
  rows.forEach((row, i) => {
    const y = 286 + i * 65;
    body += roundedRect(96, y, 245, 46, row[3], row[4], 8);
    body += textBlock(row[0], 118, y + 29, 190, { size: 15, bold: true, color: row[4], maxChars: 28 });
    body += textBlock(row[1], 392, y + 34, 150, { size: 30, bold: true, color: row[4] });
    body += textBlock(row[2], 560, y + 30, 520, { size: 14, color: C.muted, maxChars: 76 });
    if (i < rows.length - 1) {
      body += `<line x1="218" y1="${y + 48}" x2="218" y2="${y + 62}" stroke="#94A3B8" stroke-width="2"/>`;
      body += `<polygon points="218,${y + 66} 213,${y + 57} 223,${y + 57}" fill="#94A3B8"/>`;
    }
  });
  body += roundedRect(780, 608, 340, 48, C.orangeFill, "#F4C98A");
  body += textBlock("Do not compare denominators as if they are the same object.", 802, 636, 296, { size: 13, bold: true, color: C.orange, maxChars: 50 });
  return `<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="720">${body}</svg>`;
}

function slide4Svg() {
  let body = base(
    "Clinical annotation",
    "How an episode becomes UTI in this project",
    "Primary UTI status requires culture support plus catheter-aware compatible symptoms; otherwise the episode is Not_UTI.",
    4,
  );

  body += roundedRect(74, 178, 300, 160, C.blueFill, "#B9D5F0");
  body += textBlock("1. Culture support", 96, 214, 250, { size: 20, bold: true, color: C.blue });
  body += textBlock("Urine culture evidence supports possible infection under the primary lower-threshold rule.", 96, 252, 240, { size: 13.5, color: C.muted, maxChars: 38, lineHeight: 18 });
  body += textBlock("Primary rule uses >=10^3 CFU/mL support where CFU data are available.", 96, 310, 240, { size: 11.5, color: C.blue, maxChars: 42, lineHeight: 15 });

  body += roundedRect(490, 178, 300, 160, C.greenFill, "#BCE4C9");
  body += textBlock("2. Symptom rule", 512, 214, 250, { size: 20, bold: true, color: C.green });
  body += textBlock("Symptoms must be compatible with UTI, using different logic for catheter and non-catheter episodes.", 512, 252, 240, { size: 13.5, color: C.muted, maxChars: 38, lineHeight: 18 });
  body += textBlock("Non-catheter: local urinary symptoms / flank-pain context. Catheter: systemic-sign route.", 512, 310, 240, { size: 11.5, color: C.green, maxChars: 42, lineHeight: 15 });

  body += arrow(382, 258, 476, 258);
  body += textBlock("+", 425, 250, 30, { size: 26, bold: true, color: C.ink });
  body += arrow(798, 258, 890, 258);
  body += roundedRect(902, 198, 250, 118, C.orangeFill, "#F4B15F");
  body += textBlock("UTI", 1016, 244, 80, { size: 32, bold: true, color: C.orange, anchor: "middle" });
  body += textBlock("Both criteria met", 1027, 276, 160, { size: 13, color: C.orange, anchor: "middle" });

  body += roundedRect(96, 428, 500, 110, C.orangeFill, "#F4B15F");
  body += textBlock("Primary clinical denominator", 122, 462, 430, { size: 18, bold: true, color: C.orange });
  body += textBlock("583 clinical episodes = 18 UTI + 565 Not_UTI", 122, 500, 430, { size: 20, bold: true, color: C.ink, maxChars: 60 });

  body += roundedRect(682, 428, 500, 110, C.slateFill, C.line);
  body += textBlock("Not_UTI is a mixed comparator", 708, 462, 430, { size: 18, bold: true, color: C.ink });
  body += textBlock("Includes episodes without both UTI criteria: bacteriuria without symptom rule, culture-negative/below-threshold, and sensitivity-only near-miss rows.", 708, 498, 420, { size: 13, color: C.muted, maxChars: 72, lineHeight: 17 });

  body += roundedRect(202, 574, 876, 54, "#FFFFFF", C.line);
  body += textBlock("How to say it live: UTI is not just bacteria in urine; it is culture support plus a compatible clinical picture.", 232, 607, 816, { size: 15, bold: true, color: C.ink, maxChars: 112 });

  return `<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="720">${body}</svg>`;
}

async function writeSvgPng(svg, out) {
  await sharp(Buffer.from(svg)).png().toFile(out);
}

function padded(n) {
  return String(n).padStart(2, "0");
}

async function buildFinalImages() {
  await fsp.mkdir(FINAL_DIR, { recursive: true });
  const finalPaths = [];

  const slide1 = path.join(BASE_DIR, "slide-01.png");
  if (!fs.existsSync(slide1)) throw new Error(`Missing base slide: ${slide1}`);
  let out = path.join(FINAL_DIR, "slide-01.png");
  await fsp.copyFile(slide1, out);
  finalPaths.push(out);

  out = path.join(FINAL_DIR, "slide-02.png");
  await writeSvgPng(slide2Svg(), out);
  finalPaths.push(out);

  out = path.join(FINAL_DIR, "slide-03.png");
  await writeSvgPng(slide3Svg(), out);
  finalPaths.push(out);

  out = path.join(FINAL_DIR, "slide-04.png");
  await writeSvgPng(slide4Svg(), out);
  finalPaths.push(out);

  for (let old = 2; old <= 22; old += 1) {
    const src = path.join(BASE_DIR, `slide-${padded(old)}.png`);
    if (!fs.existsSync(src)) throw new Error(`Missing base slide: ${src}`);
    const dest = path.join(FINAL_DIR, `slide-${padded(old + 3)}.png`);
    await fsp.copyFile(src, dest);
    finalPaths.push(dest);
  }
  return finalPaths;
}

async function buildPptx(finalPaths) {
  await fsp.mkdir(path.dirname(OUTPUT), { recursive: true });
  const pptx = new pptxgen();
  pptx.defineLayout({ name: "CUSTOM_WIDE", width: 13.333333, height: 7.5 });
  pptx.layout = "CUSTOM_WIDE";
  pptx.author = "Aamir / Codex";
  pptx.company = "rUTI VF analysis project";
  pptx.subject = "PowerPoint-safe image deck with newcomer onboarding slides";
  pptx.title = "Longitudinal urinary E. coli VF pipeline review with onboarding";
  pptx.lang = "en-US";
  pptx.theme = { headFontFace: "Aptos Display", bodyFontFace: "Aptos", lang: "en-US" };
  for (const img of finalPaths) {
    const slide = pptx.addSlide();
    slide.background = { color: "FFFFFF" };
    slide.addImage({ path: img, x: 0, y: 0, w: 13.333333, h: 7.5 });
  }
  await pptx.writeFile({ fileName: OUTPUT });
}

async function buildContact(finalPaths) {
  const contact = spawnSync(
    PYTHON,
    [path.join(SKILL_DIR, "scripts", "make_contact_sheet.py"), "--output", CONTACT_SHEET, ...finalPaths],
    { encoding: "utf8" },
  );
  if (contact.status !== 0) throw new Error([contact.stdout, contact.stderr].filter(Boolean).join("\n"));
}

(async () => {
  const finalPaths = await buildFinalImages();
  if (finalPaths.length !== 25) throw new Error(`Expected 25 final slide images, got ${finalPaths.length}`);
  await buildPptx(finalPaths);
  await buildContact(finalPaths);
  console.log(OUTPUT);
  console.log(CONTACT_SHEET);
})().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

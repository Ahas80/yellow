#!/usr/bin/env node

import fs from "fs/promises";
import path from "path";
import { spawnSync } from "child_process";

import {
  ensureArtifactToolWorkspace,
  importArtifactTool,
  padSlideNumber,
  saveBlobToFile,
} from "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.521.10419/skills/presentations/scripts/artifact_tool_utils.mjs";

const WORKSPACE = "/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v5";
const PPTX = path.join(WORKSPACE, "output", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx");
const PREVIEW_DIR = path.join(WORKSPACE, "preview", "editable");
const LAYOUT_DIR = path.join(WORKSPACE, "layout", "editable");
const CONTACT_SHEET = path.join(WORKSPACE, "preview", "Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_contact_sheet.png");
const SKILL_DIR = "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.521.10419/skills/presentations";
const PYTHON = "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";

function slidesFromPresentation(presentation) {
  if (Array.isArray(presentation.slides?.items)) return presentation.slides.items;
  return Array.from({ length: presentation.slides.count }, (_, index) => presentation.slides.getItem(index));
}

async function main() {
  await ensureArtifactToolWorkspace(WORKSPACE);
  await fs.mkdir(PREVIEW_DIR, { recursive: true });
  await fs.mkdir(LAYOUT_DIR, { recursive: true });

  const { FileBlob, PresentationFile } = await importArtifactTool(WORKSPACE);
  const presentation = await PresentationFile.importPptx(await FileBlob.load(PPTX));
  const slides = slidesFromPresentation(presentation);
  if (slides.length !== 17) throw new Error(`Expected 17 slides, found ${slides.length}`);

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
    PYTHON,
    [path.join(SKILL_DIR, "scripts", "make_contact_sheet.py"), "--output", CONTACT_SHEET, ...previewPaths],
    { encoding: "utf8" },
  );
  if (contact.status !== 0) {
    throw new Error([contact.stdout, contact.stderr].filter(Boolean).join("\n"));
  }
  console.log(CONTACT_SHEET);
}

main().catch((error) => {
  console.error(error.stack || error.message || String(error));
  process.exit(1);
});

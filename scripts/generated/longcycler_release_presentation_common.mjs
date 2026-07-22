import crypto from "node:crypto";
import { execFile } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const execFileAsync = promisify(execFile);

export const SKILL_DIR = "/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/presentations/26.709.11516/skills/presentations";
export const PYTHON = "/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";
export const DEFAULT_PROJECT_ROOT = "/Users/Aamir/Desktop/rUTIs";

const FORBIDDEN_TOKEN = ["fl", "ye"].join("");
const TEXT_EXTENSIONS = new Set([
  ".csv", ".json", ".md", ".ndjson", ".txt", ".xml", ".rels", ".yaml", ".yml",
]);

export function parseArgs(argv) {
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

export async function run(command, args, options = {}) {
  const { stdout, stderr } = await execFileAsync(command, args, {
    cwd: options.cwd,
    env: { ...process.env, ...(options.env ?? {}) },
    encoding: options.encoding ?? "utf8",
    maxBuffer: options.maxBuffer ?? 128 * 1024 * 1024,
  });
  if (options.echo && stdout?.trim()) console.log(stdout.trim());
  if (options.echo && stderr?.trim()) console.error(stderr.trim());
  return { stdout, stderr };
}

export async function bootstrapToScratch({ scriptUrl, taskSlug, argv }) {
  if (process.env.RUTI_PRESENTATION_REEXEC === "1") return false;

  const args = parseArgs(argv);
  const projectRoot = path.resolve(args["project-root"] || DEFAULT_PROJECT_ROOT);
  const workspace = path.resolve(
    args.workspace ||
      path.join(
        os.tmpdir(),
        "codex-presentations",
        process.env.CODEX_THREAD_ID || "manual-longcycler-release",
        taskSlug,
        "tmp",
      ),
  );
  const sourceDir = path.join(workspace, "source");
  await fs.mkdir(sourceDir, { recursive: true });
  await run(process.execPath, [
    path.join(SKILL_DIR, "container_tools", "setup_artifact_tool_workspace.mjs"),
    "--workspace",
    workspace,
  ]);

  const sourceScript = fileURLToPath(scriptUrl);
  const copiedScript = path.join(sourceDir, path.basename(sourceScript));
  const commonSource = fileURLToPath(import.meta.url);
  const copiedCommon = path.join(sourceDir, path.basename(commonSource));
  await fs.copyFile(sourceScript, copiedScript);
  await fs.copyFile(commonSource, copiedCommon);

  const forwarded = [...argv];
  if (!args["project-root"]) forwarded.push("--project-root", projectRoot);
  if (!args.workspace) forwarded.push("--workspace", workspace);
  const result = await run(process.execPath, [copiedScript, ...forwarded], {
    cwd: projectRoot,
    env: { RUTI_PRESENTATION_REEXEC: "1" },
    echo: true,
  });
  if (result.stderr?.trim()) console.error(result.stderr.trim());
  return true;
}

function devRegistry() {
  return {
    schema_version: "longcycler_release_claim_registry_v1",
    generated_at: "2026-07-13 12:00:00 CEST",
    analysis_scope: {
      assembly_policy: "selected QC-passing Longcycler only",
      clinical_phenotype: "operational UTI phenotype",
      clinical_definition_version: "operational-2026-07",
      interpretation: "exploratory observational analysis; no causal claim",
    },
    method_contract: {
      operational_phenotype: {
        culture_lower_bound_cfu_per_ml: 1000,
        rule: "versioned operational culture-plus-compatible-symptom phenotype",
        caveat: "not a reconstruction of the full published protocol",
      },
      assembly_qc: {
        max_contigs: 200,
        min_n50_bp: 20000,
        min_genome_size_bp: 4000000,
        max_genome_size_bp: 6000000,
        excluded_metrics: ["read coverage", "completeness", "contamination"],
      },
      vfdb: {
        tool: "ABRicate",
        database: "VFDB",
        min_identity_pct: 80,
        min_coverage_pct: 80,
        provenance: "SHA-bound calls from the selected Longcycler FASTA manifest",
      },
      mlst: {
        role: "lineage context; not pair-specific continuity proof",
        provider_min_good_targets_pct: 95,
        provider_policy: "provider_qc95 call key/path-linked to the selected Longcycler episode; local fallback excluded",
        fallback: "labelled local MLST from the same selected Longcycler FASTA where required",
      },
      direct_pair_evidence: {
        tool: "dnadiff",
        role: "primary pair-specific distance evidence",
        operational_snp_threshold: 25,
        priority: "graph connectivity and MLST agreement cannot override a conflicting direct pair",
      },
      population_context: {
        core_genome_tool: "Parsnp",
        pangenome_tool: "Panaroo",
        role: "population context; not a substitute for direct pair evidence",
      },
    },
    analytical_cohort: {
      episodes: 532,
      residents: 161,
      operational_UTI: 16,
      operational_Not_UTI: 516,
    },
    attrition_qc_context: {
      label: "full clinical source retained only for attrition/QC context",
      episodes: 583,
      residents: 166,
      operational_UTI: 18,
      operational_Not_UTI: 565,
    },
    direct_pairs: { all_within_resident: 893 },
    adjacent_transitions: {
      pairs: 371,
      residents: 139,
      operational_snp_threshold: 25,
      at_or_below_threshold: 140,
      Not_UTI_to_UTI: 9,
      Not_UTI_to_UTI_at_or_below_threshold: 5,
    },
    mechanism_casebook: { cases: 9, linked: 9, missing: 0 },
    near_miss_audit: { rows: 17, label: "near-miss rows; not operational UTI cases" },
    selected_uti_event_genomes: { genomes: 32, residents: 29, operational_UTI: 15, operational_Not_UTI: 17 },
    genomic_dimensions: {
      VFDB_binary_features: 227,
      MLST_typed_episodes: 514,
      distinct_preferred_ST_labels: 80,
    },
    research_questions: { first: "RQ01", last: "RQ10", count: 10, retired_questions: 0 },
    plot_files: [],
    sources: [],
  };
}

export function validateRegistry(registry) {
  const a = registry?.analytical_cohort ?? {};
  const x = registry?.attrition_qc_context ?? {};
  const t = registry?.adjacent_transitions ?? {};
  const c = registry?.mechanism_casebook ?? {};
  const n = registry?.near_miss_audit ?? {};
  const r = registry?.research_questions ?? {};
  const scope = registry?.analysis_scope ?? {};
  const methods = registry?.method_contract ?? {};
  const phenotype = methods.operational_phenotype ?? {};
  const assemblyQc = methods.assembly_qc ?? {};
  const vfdb = methods.vfdb ?? {};
  const mlst = methods.mlst ?? {};
  const direct = methods.direct_pair_evidence ?? {};
  const population = methods.population_context ?? {};
  const errors = [];
  if (registry?.schema_version !== "longcycler_release_claim_registry_v1") errors.push("schema version");
  if (scope.assembly_policy !== "selected QC-passing Longcycler only") errors.push("assembly scope");
  if (scope.clinical_phenotype !== "operational UTI phenotype") errors.push("clinical phenotype");
  if (scope.interpretation !== "exploratory observational analysis; no causal claim") errors.push("interpretation boundary");
  if (phenotype.culture_lower_bound_cfu_per_ml !== 1000 || !String(phenotype.rule || "").includes("operational")) errors.push("operational phenotype method");
  if (assemblyQc.max_contigs !== 200 || assemblyQc.min_n50_bp !== 20000 || assemblyQc.min_genome_size_bp !== 4000000 || assemblyQc.max_genome_size_bp !== 6000000) errors.push("assembly QC method");
  if (vfdb.tool !== "ABRicate" || vfdb.database !== "VFDB" || vfdb.min_identity_pct !== 80 || vfdb.min_coverage_pct !== 80 || !String(vfdb.provenance || "").includes("selected Longcycler")) errors.push("VFDB method");
  if (mlst.provider_min_good_targets_pct !== 95 || !String(mlst.provider_policy || "").includes("selected Longcycler") || !String(mlst.role || "").includes("lineage context")) errors.push("MLST method");
  if (direct.tool !== "dnadiff" || direct.operational_snp_threshold !== 25 || !String(direct.priority || "").includes("cannot override")) errors.push("direct-pair method");
  if (population.core_genome_tool !== "Parsnp" || population.pangenome_tool !== "Panaroo" || !String(population.role || "").includes("population context")) errors.push("population-context methods");
  if (a.episodes !== 532 || a.residents !== 161 || a.operational_UTI !== 16 || a.operational_Not_UTI !== 516) errors.push("analytical cohort");
  if (x.episodes !== 583 || x.residents !== 166 || x.operational_UTI !== 18 || x.operational_Not_UTI !== 565 || !String(x.label || "").includes("attrition/QC context")) errors.push("attrition/QC context");
  if (registry?.direct_pairs?.all_within_resident !== 893) errors.push("direct pairs");
  if (t.pairs !== 371 || t.residents !== 139 || t.operational_snp_threshold !== 25 || t.at_or_below_threshold !== 140) errors.push("adjacent transitions");
  if (t.Not_UTI_to_UTI !== 9 || t.Not_UTI_to_UTI_at_or_below_threshold !== 5) errors.push("focused transitions");
  if (c.cases !== 9 || c.linked !== 9 || c.missing !== 0) errors.push("mechanism casebook");
  if (n.rows !== 17) errors.push("near-miss audit");
  if (r.first !== "RQ01" || r.last !== "RQ10" || r.count !== 10 || r.retired_questions !== 0) errors.push("research-question layer");
  const serialized = JSON.stringify(registry);
  if (serialized.toLowerCase().includes(FORBIDDEN_TOKEN)) errors.push("retired input token");
  if (errors.length) throw new Error(`Claim registry contract failed: ${errors.join(", ")}`);
}

export async function loadRegistry({ projectRoot, registryPath, devFixture }) {
  const resolved = registryPath
    ? path.resolve(registryPath)
    : path.join(projectRoot, "results", "pipeline", "longcycler_release_claim_registry.json");
  const registryBytes = devFixture
    ? Buffer.from(`${JSON.stringify(devRegistry(), null, 2)}\n`, "utf8")
    : await fs.readFile(resolved);
  const registry = JSON.parse(registryBytes.toString("utf8"));
  validateRegistry(registry);
  return {
    registry,
    registryPath: devFixture ? "development fixture generated in external scratch" : resolved,
    // Hash the exact bytes that were parsed. Re-serialising JSON can change
    // whitespace and therefore does not prove which registry file was used.
    registrySha256: crypto.createHash("sha256").update(registryBytes).digest("hex"),
  };
}

export function outputDirectory(projectRoot, outputRoot, family) {
  return outputRoot
    ? path.join(path.resolve(outputRoot), family)
    : path.join(projectRoot, "outputs", family);
}

export function releaseDate(registry) {
  const raw = String(registry.generated_at || "");
  const match = raw.match(/(20\d{2})-(\d{2})-(\d{2})/);
  if (!match) return "release date recorded in registry";
  const date = new Date(`${match[1]}-${match[2]}-${match[3]}T12:00:00Z`);
  return new Intl.DateTimeFormat("en-GB", {
    day: "numeric",
    month: "long",
    year: "numeric",
    timeZone: "UTC",
  }).format(date);
}

export function sha256Text(text) {
  return crypto.createHash("sha256").update(text).digest("hex");
}

export async function sha256File(filePath) {
  const bytes = await fs.readFile(filePath);
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

export async function writeTextAtomic(filePath, value) {
  await fs.mkdir(path.dirname(filePath), { recursive: true });
  const temp = path.join(path.dirname(filePath), `.${path.basename(filePath)}.${process.pid}.tmp`);
  await fs.writeFile(temp, value, "utf8");
  await fs.rename(temp, filePath);
}

export async function writeJsonAtomic(filePath, value) {
  await writeTextAtomic(filePath, `${JSON.stringify(value, null, 2)}\n`);
}

export function csvEscape(value) {
  const text = String(value ?? "");
  return /[",\n\r]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export async function writeCsvAtomic(filePath, rows, columns) {
  const lines = [columns.join(",")];
  for (const row of rows) lines.push(columns.map((column) => csvEscape(row[column])).join(","));
  await writeTextAtomic(filePath, `${lines.join("\n")}\n`);
}

export async function copyFileAtomic(source, target) {
  await fs.mkdir(path.dirname(target), { recursive: true });
  const temp = path.join(path.dirname(target), `.${path.basename(target)}.${process.pid}.tmp`);
  await fs.copyFile(source, temp);
  await fs.rename(temp, target);
}

export async function prepareTemplate({
  artifact,
  projectRoot,
  workspace,
  templatePptx,
  deckKind,
  expectedSlides,
  registryPath,
}) {
  await fs.rm(workspace, { recursive: true, force: true });
  await fs.mkdir(workspace, { recursive: true });
  const inspectScript = path.join(SKILL_DIR, "template_following_scripts", "inspect_template_deck.mjs");
  const validateScript = path.join(SKILL_DIR, "template_following_scripts", "validate_template_plan.mjs");
  const prepareScript = path.join(SKILL_DIR, "template_following_scripts", "prepare_template_starter_deck.mjs");
  const mapScript = path.join(projectRoot, "scripts", "generated", "create_lecturer_template_map.mjs");
  await run(process.execPath, [inspectScript, "--workspace", workspace, "--pptx", templatePptx, "--scale", "1.5"]);
  const sourcePresentation = await artifact.PresentationFile.importPptx(await artifact.FileBlob.load(templatePptx));
  const completeInspect = await sourcePresentation.inspect({
    kind: "slide,textbox,shape,image,table,chart",
    maxChars: 2000000,
  });
  await fs.writeFile(
    path.join(workspace, "template-inspect", "template-inspect.ndjson"),
    completeInspect.ndjson ?? "",
    "utf8",
  );
  const manifestPath = path.join(workspace, "template-inspect", "template-manifest.json");
  const manifest = JSON.parse(await fs.readFile(manifestPath, "utf8"));
  manifest.inspectTruncated = Boolean(completeInspect.truncated);
  manifest.inspectMetadata = completeInspect.metadata ?? {};
  await fs.writeFile(manifestPath, `${JSON.stringify(manifest, null, 2)}\n`, "utf8");
  await run(process.execPath, [
    mapScript,
    "--workspace", workspace,
    "--deck-kind", deckKind,
    "--expected-slides", String(expectedSlides),
    "--registry", registryPath,
  ]);
  await run(process.execPath, [validateScript, "--workspace", workspace, "--map", path.join(workspace, "template-frame-map.json")]);
  await run(process.execPath, [
    prepareScript,
    "--workspace", workspace,
    "--pptx", templatePptx,
    "--map", path.join(workspace, "template-frame-map.json"),
    "--out", path.join(workspace, "template-starter.pptx"),
    "--preview-dir", path.join(workspace, "template-starter-preview"),
    "--layout-dir", path.join(workspace, "template-starter-layout"),
    "--contact-sheet", path.join(workspace, "template-starter-contact-sheet.png"),
    "--scale", "1.5",
  ], { env: { PYTHON } });
}

async function writeBlob(filePath, blob) {
  await fs.mkdir(path.dirname(filePath), { recursive: true });
  await fs.writeFile(filePath, new Uint8Array(await blob.arrayBuffer()));
}

function slidesFromPresentation(presentation) {
  if (Array.isArray(presentation.slides?.items)) return presentation.slides.items;
  if (Number.isInteger(presentation.slides?.count) && typeof presentation.slides.getItem === "function") {
    return Array.from({ length: presentation.slides.count }, (_, index) => presentation.slides.getItem(index));
  }
  throw new Error("Could not enumerate imported slides.");
}

export async function editTemplateDeck({
  artifact,
  projectRoot,
  workspace,
  templatePptx,
  outputPptx,
  deckKind,
  textBySlide,
  notesBySlide = {},
  registryPath,
}) {
  const expectedSlides = Object.keys(textBySlide).length;
  await prepareTemplate({ artifact, projectRoot, workspace, templatePptx, deckKind, expectedSlides, registryPath });
  await fs.copyFile(fileURLToPath(import.meta.url), path.join(workspace, "authoring-source.mjs"));

  const { FileBlob, PresentationFile } = artifact;
  const presentation = await PresentationFile.importPptx(
    await FileBlob.load(path.join(workspace, "template-starter.pptx")),
  );
  const slides = slidesFromPresentation(presentation);
  if (slides.length !== expectedSlides) throw new Error(`${deckKind}: expected ${expectedSlides} slides; imported ${slides.length}`);

  const previewDir = path.join(workspace, "final-preview");
  const layoutDir = path.join(workspace, "layout", "final");
  await fs.mkdir(previewDir, { recursive: true });
  await fs.mkdir(layoutDir, { recursive: true });

  for (let slideNumber = 1; slideNumber <= slides.length; slideNumber += 1) {
    const slide = slides[slideNumber - 1];
    const padded = String(slideNumber).padStart(2, "0");
    // Artifact IDs are scoped to the in-memory import. The IDs recorded while
    // preparing template-starter.pptx are not guaranteed to survive the next
    // import, so resolve the inherited objects from this live slide layout.
    const importedLayout = await slide.export({ format: "layout" });
    const layout = JSON.parse(await importedLayout.text());
    const textElements = layout.elements.filter((element) => typeof element.text === "string");
    const replacements = textBySlide[slideNumber];
    if (!Array.isArray(replacements) || replacements.length !== textElements.length) {
      throw new Error(`${deckKind} slide ${slideNumber}: ${textElements.length} inherited text objects but ${replacements?.length ?? 0} replacements`);
    }

    for (let index = 0; index < textElements.length; index += 1) {
      const element = textElements[index];
      if (typeof element.aid !== "string" || !element.aid) {
        throw new Error(`${deckKind} slide ${slideNumber}: inherited text element lacks an artifact ID`);
      }
      const target = presentation.resolve(element.aid);
      if (!target?.text) throw new Error(`${deckKind} slide ${slideNumber}: unresolved inherited text anchor ${element.aid}`);
      target.text.set(String(replacements[index]));
    }

    if (notesBySlide[slideNumber]) {
      slide.speakerNotes.textFrame.setText(notesBySlide[slideNumber]);
      slide.speakerNotes.setVisible(true);
    }

    await writeBlob(
      path.join(previewDir, `slide-${padded}.png`),
      await presentation.export({ slide, format: "png", scale: 2 }),
    );
    const layoutBlob = await slide.export({ format: "layout" });
    await writeBlob(path.join(layoutDir, `slide-${padded}.layout.json`), layoutBlob);
  }

  await writeBlob(
    path.join(workspace, "final-montage.webp"),
    await presentation.export({ format: "webp", montage: true, scale: 1 }),
  );
  const inspect = await presentation.inspect({
    kind: "slide,textbox,shape,notes",
    maxChars: 300000,
  });
  await fs.writeFile(path.join(workspace, "final-inspect.ndjson"), inspect.ndjson ?? "", "utf8");

  await fs.mkdir(path.dirname(outputPptx), { recursive: true });
  const tempPptx = path.join(path.dirname(outputPptx), `.${path.basename(outputPptx)}.${process.pid}.tmp.pptx`);
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(tempPptx);
  await assertNoEmptyPlaceholders(tempPptx);
  await assertPptxHasNoForbiddenContent(tempPptx);
  await fs.rename(tempPptx, outputPptx);

  await run(process.execPath, [
    path.join(SKILL_DIR, "template_following_scripts", "check_template_fidelity.mjs"),
    "--workspace", workspace,
    "--starter-pptx", path.join(workspace, "template-starter.pptx"),
    "--final-pptx", outputPptx,
    "--map", path.join(workspace, "template-frame-map.json"),
    "--starter-layout-dir", path.join(workspace, "template-starter-layout"),
    "--final-layout-dir", layoutDir,
    "--edit-dir", workspace,
  ]);
  await run(PYTHON, [path.join(SKILL_DIR, "container_tools", "slides_test.py"), outputPptx]);
  return {
    outputPptx,
    workspace,
    previewDir,
    layoutDir,
    inspectPath: path.join(workspace, "final-inspect.ndjson"),
    montagePath: path.join(workspace, "final-montage.webp"),
    slideCount: slides.length,
  };
}

async function zipEntryNames(pptxPath) {
  const { stdout } = await run("unzip", ["-Z1", pptxPath], { maxBuffer: 32 * 1024 * 1024 });
  return stdout.split(/\r?\n/).filter(Boolean);
}

async function readZipEntry(pptxPath, entry) {
  const unzipPattern = entry.replace(/[\\[\]*?]/g, "\\$&");
  const { stdout } = await run("unzip", ["-p", pptxPath, unzipPattern], {
    encoding: "buffer",
    maxBuffer: 128 * 1024 * 1024,
  });
  return Buffer.from(stdout).toString("utf8");
}

export async function assertNoEmptyPlaceholders(pptxPath) {
  const entries = (await zipEntryNames(pptxPath)).filter((entry) => /^ppt\/slides\/slide\d+\.xml$/.test(entry));
  const failures = [];
  for (const entry of entries) {
    const xml = await readZipEntry(pptxPath, entry);
    for (const shape of xml.matchAll(/<p:sp\b[\s\S]*?<\/p:sp>/g)) {
      const body = shape[0];
      if (!/<p:ph\b/.test(body)) continue;
      const text = [...body.matchAll(/<a:t>([\s\S]*?)<\/a:t>/g)]
        .map((match) => match[1].replace(/<[^>]+>/g, "").trim())
        .join("")
        .trim();
      if (!text) failures.push(entry);
    }
  }
  if (failures.length) throw new Error(`Empty inherited placeholders remain in ${pptxPath}: ${[...new Set(failures)].join(", ")}`);
}

export async function assertPptxHasNoForbiddenContent(pptxPath) {
  const entries = (await zipEntryNames(pptxPath)).filter((entry) => /\.(?:xml|rels|txt|json)$/i.test(entry));
  for (const entry of entries) {
    const text = await readZipEntry(pptxPath, entry);
    if (text.toLowerCase().includes(FORBIDDEN_TOKEN)) {
      throw new Error(`Retired input token remains in ${path.basename(pptxPath)} package part ${entry}`);
    }
  }
}

async function listFiles(root) {
  const out = [];
  for (const entry of await fs.readdir(root, { withFileTypes: true })) {
    const full = path.join(root, entry.name);
    if (entry.isDirectory()) out.push(...(await listFiles(full)));
    else if (entry.isFile()) out.push(full);
  }
  return out;
}

export async function assertOutputTreeClean(root) {
  const files = await listFiles(root);
  for (const file of files) {
    if (file.toLowerCase().includes(FORBIDDEN_TOKEN)) throw new Error(`Retired input token remains in output path: ${file}`);
    if (path.extname(file).toLowerCase() === ".pptx") {
      await assertPptxHasNoForbiddenContent(file);
      continue;
    }
    if (TEXT_EXTENSIONS.has(path.extname(file).toLowerCase())) {
      const text = await fs.readFile(file, "utf8");
      if (text.toLowerCase().includes(FORBIDDEN_TOKEN)) throw new Error(`Retired input token remains in generated file: ${file}`);
    }
  }
}

export function releaseCounts(registry) {
  const a = registry.analytical_cohort;
  const x = registry.attrition_qc_context;
  const t = registry.adjacent_transitions;
  const g = registry.genomic_dimensions;
  const c = registry.mechanism_casebook;
  const methods = registry.method_contract;
  return {
    episodes: a.episodes,
    residents: a.residents,
    uti: a.operational_UTI,
    notUti: a.operational_Not_UTI,
    sourceEpisodes: x.episodes,
    sourceResidents: x.residents,
    sourceUti: x.operational_UTI,
    sourceNotUti: x.operational_Not_UTI,
    directPairs: registry.direct_pairs.all_within_resident,
    transitions: t.pairs,
    transitionResidents: t.residents,
    snpThreshold: t.operational_snp_threshold,
    transitionsAtThreshold: t.at_or_below_threshold,
    toUti: t.Not_UTI_to_UTI,
    toUtiAtThreshold: t.Not_UTI_to_UTI_at_or_below_threshold,
    casebookCases: c.cases,
    casebookLinked: c.linked,
    casebookMissing: c.missing,
    nearMiss: registry.near_miss_audit.rows,
    vfFeatures: g.VFDB_binary_features,
    mlstTyped: g.MLST_typed_episodes,
    mlstDistinct: g.distinct_preferred_ST_labels,
    cultureLowerBoundCfu: methods.operational_phenotype.culture_lower_bound_cfu_per_ml,
    assemblyMaxContigs: methods.assembly_qc.max_contigs,
    assemblyMinN50: methods.assembly_qc.min_n50_bp,
    assemblyMinGenomeSize: methods.assembly_qc.min_genome_size_bp,
    assemblyMaxGenomeSize: methods.assembly_qc.max_genome_size_bp,
    vfdbTool: methods.vfdb.tool,
    vfdbDatabase: methods.vfdb.database,
    vfdbMinIdentity: methods.vfdb.min_identity_pct,
    vfdbMinCoverage: methods.vfdb.min_coverage_pct,
    providerMlstMinGoodTargets: methods.mlst.provider_min_good_targets_pct,
    providerMlstPolicy: methods.mlst.provider_policy,
    directPairTool: methods.direct_pair_evidence.tool,
    coreGenomeTool: methods.population_context.core_genome_tool,
    pangenomeTool: methods.population_context.pangenome_tool,
  };
}

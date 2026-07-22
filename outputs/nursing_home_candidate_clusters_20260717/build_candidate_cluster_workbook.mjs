import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

process.on("uncaughtException", (error) => {
  console.error(`WORKBOOK_BUILD_ERROR: ${error?.message ?? error}`);
  process.exit(1);
});
process.on("unhandledRejection", (error) => {
  console.error(`WORKBOOK_BUILD_ERROR: ${error?.message ?? error}`);
  process.exit(1);
});

const baseDir = path.resolve("outputs/nursing_home_candidate_clusters_20260717");
const outputPath = path.join(baseDir, "Genetics_first_candidate_clusters.xlsx");
const previewDir = path.join(baseDir, "qa_workbook");
await fs.mkdir(previewDir, { recursive: true });

const summary = JSON.parse(
  await fs.readFile(path.join(baseDir, "analysis_summary.json"), "utf8"),
);

function parseCSV(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i += 1;
      } else if (ch === '"') {
        quoted = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows.filter((r) => r.some((v) => v !== ""));
}

function colLetter(indexZeroBased) {
  let n = indexZeroBased + 1;
  let out = "";
  while (n > 0) {
    const rem = (n - 1) % 26;
    out = String.fromCharCode(65 + rem) + out;
    n = Math.floor((n - 1) / 26);
  }
  return out;
}

function isBooleanHeader(header) {
  return /(^QC_PASS$|reliable|identical|detected|available|flag$|any_selected|any_QC|include_primary|expected_include)/i.test(
    header,
  );
}

function isDateHeader(header) {
  return /(^|_)date(_|$)|earliest_collection|latest_collection/i.test(header);
}

function isNumericHeader(header) {
  if (/(^|_)(id|ids|path|hash|timepoint|episode|sample_key)($|_)/i.test(header)) {
    return false;
  }
  if (/^ST$|^ST_[AB]$|^Participant/i.test(header)) return false;
  return /(^n_|_n$|_count|_days$|span_days|gap_days|snps|identity|jaccard|PercGoodTargets|intersection|union|min$|max$)/i.test(
    header,
  );
}

function typedValue(header, raw) {
  if (raw === "") return null;
  if (isBooleanHeader(header) && /^(TRUE|FALSE)$/i.test(raw)) {
    return raw.toUpperCase() === "TRUE";
  }
  if (isDateHeader(header) && /^\d{4}-\d{2}-\d{2}$/.test(raw)) {
    return new Date(`${raw}T00:00:00Z`);
  }
  if (isNumericHeader(header)) {
    const n = Number(raw);
    if (Number.isFinite(n)) return n;
  }
  return raw;
}

async function loadCSV(filename) {
  const rows = parseCSV(await fs.readFile(path.join(baseDir, filename), "utf8"));
  const headers = rows[0];
  const data = rows.slice(1).map((row) =>
    headers.map((header, idx) => typedValue(header, row[idx] ?? "")),
  );
  return { headers, data };
}

const sheetSpecs = [
  {
    name: "Cluster_summary",
    file: "cluster_summary.csv",
    title: "High-priority candidate genomic clusters",
    note: `${summary.high_priority_cluster_count} clusters, ${summary.high_priority_cluster_isolate_count} isolates. Uniform high-priority criteria are stated once here; detailed preliminary SNP, VF and pair-level AMR evidence remains in Candidate_pairs. These genetics-first candidates require nursing-home, ward and temporal confirmation.`,
    freezeCols: 2,
    omitHeaders: [
      "interpretation",
      "possible_cross_participant_pair_n",
      "chain_linked_flag",
      "preliminary_core_snp_min",
      "preliminary_core_snp_max",
      "vf_jaccard_min",
      "vf_jaccard_max",
      "variable_amr_genes_excluding_mdfA",
      "mdfA_detected_in_all_members",
      "shared_amr_classes_all_members",
      "reporting_note",
    ],
    fillNoneHeaders: [
      "shared_amr_genes_all_members_excluding_mdfA",
      "resistance_classes_present",
    ],
    headerRenames: {
      high_priority_pair_n: "cross_participant_pair_count",
      shared_amr_genes_all_members_excluding_mdfA:
        "shared_amr_genes_excl_mdfA",
      resistance_classes_present:
        "shared_resistance_classes_excl_mdfA",
    },
  },
  {
    name: "Candidate_pairs",
    file: "candidate_pairs.csv",
    title: "Cross-participant candidate pairs",
    note: `${summary.preliminary_candidate_pairs} screened pairs: ${summary.priority_counts["High-priority candidate"]} high priority, ${summary.priority_counts["Moderate-priority candidate"]} moderate priority, ${summary.priority_counts["Possible related lineage"]} possible related lineage, and ${summary.priority_counts["Same lineage, not the same strain"]} same-lineage pairs above the direct 25-SNP threshold.`,
    freezeCols: 4,
  },
  {
    name: "Isolate_profiles",
    file: "isolate_profiles.csv",
    title: "QC-selected isolate profiles",
    note: `${summary.active_isolates} QC-passing isolates from ${summary.active_participants} participants. AMR findings are genotypic ResFinder calls; phenotypic susceptibility was not available.`,
    freezeCols: 4,
  },
  {
    name: "Excluded_isolates",
    file: "excluded_isolates.csv",
    title: "Rows excluded from the current genomic cohort",
    note: `${summary.excluded_emailed_rows} rows in the emailed 556-row file are outside the current 532-isolate QC-selected cohort. The reason is documented for each row.`,
    freezeCols: 3,
  },
  {
    name: "Interpretation_guide",
    file: "interpretation_guide.csv",
    title: "Interpretation and reporting guide",
    note: "Use the recommended wording to keep genomic relatedness, epidemiological linkage and predicted resistance clearly separated.",
    freezeCols: 2,
  },
];

const workbook = Workbook.create();
console.log("Workbook created");

const palette = {
  navy: "#16324F",
  blue: "#2E74B5",
  paleBlue: "#E8F1F8",
  paleGreen: "#E2F0D9",
  green: "#2F6B3A",
  paleAmber: "#FFF2CC",
  amber: "#8A5B00",
  paleRed: "#FCE4D6",
  red: "#9C0006",
  gray: "#5B6573",
  paleGray: "#F2F4F7",
  border: "#D7DEE8",
  white: "#FFFFFF",
};

function widthForHeader(header) {
  if (/^candidate_cluster_id$|^screen_group_id$|^priority_category$/i.test(header)) return 22;
  if (/Participant|Episode_ID|Isolate_ID|sample_key/i.test(header)) return 18;
  if (/Collection_Date|collection_date|earliest|latest/i.test(header)) return 14;
  if (/ST_source|validation_status/i.test(header)) return 26;
  if (/^ST($|_)/i.test(header)) return 10;
  if (/reason|interpretation|wording|note|meaning/i.test(header)) return 48;
  if (/amr|resistance|vf_genes/i.test(header)) return 38;
  if (/path/i.test(header)) return 32;
  if (/hash/i.test(header)) return 18;
  if (/jaccard|identity|snps|count|_n$|days/i.test(header)) return 14;
  if (/Event_type|Status/i.test(header)) return 18;
  return Math.min(24, Math.max(11, header.length + 2));
}

for (const spec of sheetSpecs) {
  console.log(`Building ${spec.name}`);
  const loaded = await loadCSV(spec.file);
  let headers = loaded.headers;
  let data = loaded.data;
  for (const header of spec.fillNoneHeaders ?? []) {
    const idx = headers.indexOf(header);
    if (idx >= 0) {
      data = data.map((row) =>
        row.map((value, colIdx) =>
          colIdx === idx && (value === null || value === "") ? "None detected" : value,
        ),
      );
    }
  }
  if (spec.omitHeaders?.length) {
    const keepIndices = headers
      .map((header, idx) => ({ header, idx }))
      .filter(({ header }) => !spec.omitHeaders.includes(header))
      .map(({ idx }) => idx);
    headers = keepIndices.map((idx) => headers[idx]);
    data = data.map((row) => keepIndices.map((idx) => row[idx]));
  }
  headers = headers.map((header) => spec.headerRenames?.[header] ?? header);
  const sheet = workbook.worksheets.add(spec.name);
  sheet.showGridLines = false;
  const lastCol = colLetter(headers.length - 1);
  const titleBandLastCol = colLetter(Math.min(headers.length, 8) - 1);
  const tableStartRow = 5;
  const dataStartRow = tableStartRow + 1;
  const tableEndRow = tableStartRow + data.length;

  sheet.getRange(`A1:${titleBandLastCol}1`).merge();
  sheet.getRange("A1").values = [[spec.title]];
  sheet.getRange(`A2:${titleBandLastCol}2`).merge();
  sheet.getRange("A2").values = [[spec.note]];
  sheet.getRange(`A1:${titleBandLastCol}1`).format = {
    fill: palette.navy,
    font: { name: "Aptos Display", size: 18, bold: true, color: palette.white },
    verticalAlignment: "center",
  };
  sheet.getRange(`A1:${titleBandLastCol}1`).format.rowHeight = 30;
  sheet.getRange(`A2:${titleBandLastCol}2`).format = {
    fill: palette.paleBlue,
    font: { name: "Aptos", size: 10, italic: true, color: palette.navy },
    wrapText: true,
    verticalAlignment: "center",
  };
  sheet.getRange(`A2:${titleBandLastCol}2`).format.rowHeight = 34;

  sheet.getRange(`A${tableStartRow}:${lastCol}${tableEndRow}`).values = [
    headers,
    ...data,
  ];
  const table = sheet.tables.add(
    `A${tableStartRow}:${lastCol}${tableEndRow}`,
    true,
    `${spec.name.replace(/[^A-Za-z0-9]/g, "")}Table`,
  );
  table.style = "TableStyleMedium2";

  sheet.getRange(`A${tableStartRow}:${lastCol}${tableStartRow}`).format = {
    fill: palette.blue,
    font: { name: "Aptos", size: 10, bold: true, color: palette.white },
    wrapText: true,
    verticalAlignment: "center",
    borders: { preset: "outside", style: "thin", color: palette.navy },
  };
  sheet.getRange(`A${tableStartRow}:${lastCol}${tableStartRow}`).format.rowHeight = 34;
  if (data.length) {
    sheet.getRange(`A${dataStartRow}:${lastCol}${tableEndRow}`).format = {
      font: { name: "Aptos", size: 9, color: "#1F2937" },
      verticalAlignment: "center",
      borders: {
        insideHorizontal: { style: "thin", color: palette.border },
      },
    };
  }

  headers.forEach((header, idx) => {
    const col = colLetter(idx);
    sheet.getRange(`${col}${tableStartRow}:${col}${tableEndRow}`).format.columnWidth =
      widthForHeader(header);
    if (isDateHeader(header) && data.length) {
      sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).format.numberFormat =
        "yyyy-mm-dd";
    } else if (/jaccard/i.test(header) && data.length) {
      sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).format.numberFormat =
        "0.0%";
    } else if (/identity/i.test(header) && data.length) {
      sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).format.numberFormat =
        "0.0000";
    } else if (/snps|_n$|^n_|count|days|intersection|union/i.test(header) && data.length) {
      sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).format.numberFormat =
        "#,##0";
    }
    if (/reason|interpretation|wording|note|meaning|amr|resistance|vf_genes|path/i.test(header)) {
      sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).format.wrapText = true;
    }
  });

  const priorityIdx = headers.indexOf("priority_category");
  if (priorityIdx >= 0 && data.length) {
    const col = colLetter(priorityIdx);
    const priorityRange = sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`);
    priorityRange.conditionalFormats.add("containsText", {
      text: "High-priority",
      format: { fill: palette.paleGreen, font: { bold: true, color: palette.green } },
    });
    priorityRange.conditionalFormats.add("containsText", {
      text: "Moderate-priority",
      format: { fill: palette.paleBlue, font: { bold: true, color: palette.blue } },
    });
    priorityRange.conditionalFormats.add("containsText", {
      text: "Possible related",
      format: { fill: palette.paleAmber, font: { color: palette.amber } },
    });
    priorityRange.conditionalFormats.add("containsText", {
      text: "not the same strain",
      format: { fill: palette.paleGray, font: { color: palette.gray } },
    });
  }

  const clusterIdx = headers.indexOf("candidate_cluster_id");
  if (clusterIdx >= 0 && data.length) {
    const col = colLetter(clusterIdx);
    sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).conditionalFormats.add(
      "notContainsBlanks",
      { format: { fill: palette.paleGreen, font: { bold: true, color: palette.green } } },
    );
  }

  const chainIdx = headers.indexOf("chain_linked_flag");
  if (chainIdx >= 0 && data.length) {
    const col = colLetter(chainIdx);
    sheet.getRange(`${col}${dataStartRow}:${col}${tableEndRow}`).conditionalFormats.add(
      "cellIs",
      {
        operator: "equal",
        formula: "TRUE",
        format: { fill: palette.paleRed, font: { bold: true, color: palette.red } },
      },
    );
  }

  sheet.freezePanes.freezeRows(tableStartRow);
  if (spec.freezeCols) sheet.freezePanes.freezeColumns(spec.freezeCols);
}

const clusterCheck = await workbook.inspect({
  kind: "table",
  range: "Cluster_summary!A1:O16",
  include: "values,formulas",
  tableMaxRows: 16,
  tableMaxCols: 15,
  maxChars: 8000,
});
console.log(clusterCheck.ndjson);

const errorCheck = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
  maxChars: 3000,
});
console.log(errorCheck.ndjson);

const previewRanges = {
  Cluster_summary: "A1:O16",
  Candidate_pairs: "A1:N24",
  Isolate_profiles: "A1:N24",
  Excluded_isolates: "A1:N29",
  Interpretation_guide: "A1:E24",
};
for (const [sheetName, range] of Object.entries(previewRanges)) {
  const preview = await workbook.render({
    sheetName,
    range,
    scale: 1.2,
    format: "png",
  });
  await fs.writeFile(
    path.join(previewDir, `${sheetName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`Saved ${outputPath}`);

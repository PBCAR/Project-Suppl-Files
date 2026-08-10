import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outputDir = process.argv[2];
if (!outputDir) {
  throw new Error("Usage: node build_publication_tables.mjs <v20-run-output-dir> [banner]");
}
const banner = process.argv[3] ?? "";
const tableDir = path.join(outputDir, "tables");
const previewDir = path.join(tableDir, "previews");
await fs.mkdir(previewDir, { recursive: true });

const evidencePath = path.join(tableDir, "table_1_polygenic_liability_evidence.csv");
const samplePath = path.join(
  tableDir,
  "supplementary_table_s1_analytic_sample_characteristics.csv",
);
const supplementPath = path.join(
  tableDir,
  "supplementary_table_s2_all_cannabis_phenotypes.csv",
);
const notesPath = path.join(tableDir, "baseline_table_notes.csv");

async function csvValues(csvPath) {
  const csvText = await fs.readFile(csvPath, "utf8");
  const csvWorkbook = await Workbook.fromCSV(csvText, { sheetName: "Data" });
  return csvWorkbook.worksheets.getItem("Data").getUsedRange(true).values;
}

function selectColumns(values, columns, labels) {
  const sourceHeaders = values[0].map((value) => String(value));
  const indexes = columns.map((column) => {
    const index = sourceHeaders.indexOf(column);
    if (index < 0) throw new Error(`Required CSV column missing: ${column}`);
    return index;
  });
  return [
    labels,
    ...values.slice(1).map((row) => indexes.map((index) => row[index] ?? null)),
  ];
}

const evidenceValues = selectColumns(
  await csvValues(evidencePath),
  [
    "polygenic_score",
    "pgs_domain",
    "baseline_fdr_phenotypes_n",
    "shared_with_externalizing_n",
    "outside_externalizing_n",
    "selected_for_cud_followup",
    "pgs_by_wave_p_value",
    "pgs_by_wave_fdr_q",
    "later_cud_fdr_waves",
    "cud_followup_selection_reason",
  ],
  [
    "Polygenic score",
    "Domain",
    "Baseline FDR phenotypes",
    "Shared with externalizing",
    "Outside externalizing",
    "Selected for CUD follow-up",
    "PGS-by-wave P",
    "PGS-by-wave FDR q",
    "FDR-significant CUD waves",
    "Selection reason",
  ],
);

const sampleValues = selectColumns(
  await csvValues(samplePath),
  [
    "section",
    "characteristic",
    "scale_or_coding",
    "denominator_population",
    "overall_statistic",
    "overall_n",
    "female_statistic",
    "female_n",
    "male_statistic",
    "male_n",
    "notes",
  ],
  [
    "Section",
    "Characteristic",
    "Scale / coding",
    "Population",
    "Overall",
    "Overall n",
    "Female",
    "Female n",
    "Male",
    "Male n",
    "Notes",
  ],
);

const supplementValues = selectColumns(
  await csvValues(supplementPath),
  [
    "domain",
    "characteristic",
    "family",
    "scale_or_coding",
    "denominator_population",
    "overall_statistic",
    "overall_n",
    "female_statistic",
    "female_n",
    "male_statistic",
    "male_n",
    "sex_difference_test",
    "sex_difference_n",
    "sex_difference_p",
    "sex_difference_fdr_q",
  ],
  [
    "Domain",
    "Cannabis phenotype",
    "Type",
    "Scale / coding",
    "Population",
    "Overall",
    "Overall n",
    "Female",
    "Female n",
    "Male",
    "Male n",
    "Sex comparison",
    "Test n",
    "P",
    "BH q",
  ],
);
for (const row of supplementValues.slice(1)) {
  for (const columnIndex of [13, 14]) {
    const probability = Number(row[columnIndex]);
    if (Number.isFinite(probability)) {
      row[columnIndex] = probability < 0.001 ? "<0.001" : probability.toFixed(3);
    }
  }
}

const phenotypeNotesValues = selectColumns(
  await csvValues(supplementPath),
  [
    "domain",
    "characteristic",
    "scale_or_coding",
    "denominator_population",
    "notes",
  ],
  ["Domain", "Cannabis phenotype", "Scale / coding", "Population", "Notes"],
);

const notesSource = await csvValues(notesPath);
const notesValues = selectColumns(notesSource, ["note_order", "note"], ["No.", "Table note"]);

const workbook = Workbook.create();
const teal = "#145A64";
const tealDark = "#0E3F46";
const paleA = "#EAF4F5";
const paleB = "#F6F9F9";
const border = "#C9D8DA";

function stylePublicationSheet(sheet, values, title, subtitle, widths) {
  const rowCount = values.length + 3;
  const columnCount = values[0].length;
  const lastColumn = String.fromCharCode(64 + columnCount);
  sheet.showGridLines = false;
  sheet.getRange(`A1:${lastColumn}1`).merge();
  sheet.getRange("A1").values = [[title]];
  sheet.getRange(`A2:${lastColumn}2`).merge();
  sheet.getRange("A2").values = [[banner ? `${banner} — ${subtitle}` : subtitle]];
  sheet.getRange(`A4:${lastColumn}${rowCount}`).values = values;

  sheet.getRange(`A1:${lastColumn}1`).format = {
    fill: tealDark,
    font: { bold: true, color: "#FFFFFF", size: 14 },
    verticalAlignment: "center",
    horizontalAlignment: "left",
  };
  sheet.getRange(`A1:${lastColumn}1`).format.rowHeight = 28;
  sheet.getRange(`A2:${lastColumn}2`).format = {
    fill: "#DDECEE",
    font: { italic: true, color: "#24464B", size: 10 },
    verticalAlignment: "center",
    horizontalAlignment: "left",
    wrapText: true,
  };
  sheet.getRange(`A2:${lastColumn}2`).format.rowHeight = 34;

  const header = sheet.getRange(`A4:${lastColumn}4`);
  header.format = {
    fill: teal,
    font: { bold: true, color: "#FFFFFF", size: 10 },
    verticalAlignment: "center",
    horizontalAlignment: "left",
    wrapText: true,
    borders: { preset: "outside", style: "medium", color: tealDark },
  };
  header.format.rowHeight = 34;

  const body = sheet.getRange(`A5:${lastColumn}${rowCount}`);
  body.format = {
    font: { color: "#202124", size: 10 },
    verticalAlignment: "top",
    horizontalAlignment: "left",
    wrapText: true,
    borders: {
      insideHorizontal: { style: "thin", color: border },
      bottom: { style: "thin", color: "#9FB5B9" },
    },
  };

  let currentSection = null;
  let sectionIndex = -1;
  for (let index = 1; index < values.length; index += 1) {
    const section = values[index][0];
    if (section !== currentSection) {
      currentSection = section;
      sectionIndex += 1;
    }
    sheet.getRange(`A${index + 4}:${lastColumn}${index + 4}`).format.fill =
      sectionIndex % 2 === 0 ? paleA : paleB;
  }
  sheet.getRange(`A5:B${rowCount}`).format.font = { bold: true, color: "#164E57" };
  for (let index = 0; index < widths.length; index += 1) {
    sheet.getRangeByIndexes(0, index, rowCount, 1).format.columnWidthPx = widths[index];
  }
  body.format.autofitRows();
  sheet.freezePanes.freezeRows(4);
  sheet.freezePanes.freezeColumns(2);
}

const evidenceSheet = workbook.worksheets.add("Table 1");
stylePublicationSheet(
  evidenceSheet,
  evidenceValues,
  "Table 1. Polygenic association coverage and repeated CUD evidence",
  "All 12 baseline PGSs are shown; CUD follow-up is restricted by the audited selection rule.",
  [210, 260, 120, 130, 130, 135, 105, 115, 170, 420],
);

const sampleSheet = workbook.worksheets.add("Table S1");
stylePublicationSheet(
  sampleSheet,
  sampleValues,
  "Table S1. Baseline characteristics of the analytic sample",
  "Combined and sex-stratified summaries with variable-specific denominators; no descriptive hypothesis tests.",
  [150, 230, 300, 200, 125, 85, 125, 85, 125, 85, 360],
);

const supplementSheet = workbook.worksheets.add("Table S2");
stylePublicationSheet(
  supplementSheet,
  supplementValues,
  "Table S2. Baseline distributions of all cannabis-related phenotypes",
  "Combined and sex-stratified summaries and two-sample sex comparisons use each phenotype's designated analysis population; q values apply BH correction across 30 comparisons.",
  [170, 230, 100, 340, 200, 125, 85, 125, 85, 125, 85, 190, 85, 85, 85],
);

const phenotypeNotesSheet = workbook.worksheets.add("Phenotype notes");
stylePublicationSheet(
  phenotypeNotesSheet,
  phenotypeNotesValues,
  "Supplementary phenotype definitions and scoring notes",
  "Purpose-built legalization-context perception items are identified and their scaling is reported explicitly.",
  [170, 230, 420, 220, 600],
);

const notesSheet = workbook.worksheets.add("Table notes");
notesSheet.showGridLines = false;
notesSheet.getRange("A1:B1").merge();
notesSheet.getRange("A1").values = [["Table construction and interpretation notes"]];
notesSheet.getRange(`A3:B${notesValues.length + 2}`).values = notesValues;
notesSheet.getRange("A1:B1").format = {
  fill: tealDark,
  font: { bold: true, color: "#FFFFFF", size: 14 },
};
notesSheet.getRange("A3:B3").format = {
  fill: teal,
  font: { bold: true, color: "#FFFFFF", size: 10 },
};
notesSheet.getRange(`A4:B${notesValues.length + 2}`).format = {
  fill: paleA,
  font: { color: "#202124", size: 10 },
  verticalAlignment: "top",
  wrapText: true,
  borders: { insideHorizontal: { style: "thin", color: border } },
};
notesSheet.getRange(`A1:A${notesValues.length + 2}`).format.columnWidthPx = 65;
notesSheet.getRange(`B1:B${notesValues.length + 2}`).format.columnWidthPx = 900;
notesSheet.getRange(`A3:B${notesValues.length + 2}`).format.autofitRows();
notesSheet.freezePanes.freezeRows(3);

const inspect = await workbook.inspect({
  kind: "table",
  range: "Table 1!A1:J16",
  include: "values,formulas",
  tableMaxRows: 12,
  tableMaxCols: 10,
  maxChars: 12000,
});
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 200 },
  summary: "final formula error scan",
});
await fs.writeFile(
  path.join(tableDir, "PATH_CANN_publication_tables_v20.inspect.ndjson"),
  `${inspect.ndjson}\n${errors.ndjson}\n`,
  "utf8",
);

for (const [sheetName, fileName] of [
  ["Table 1", "table_1_polygenic_liability_evidence.png"],
  ["Table S1", "supplementary_table_s1_analytic_sample_characteristics.png"],
  ["Table S2", "supplementary_table_s2_all_cannabis_phenotypes.png"],
  ["Phenotype notes", "supplementary_phenotype_notes.png"],
  ["Table notes", "table_notes.png"],
]) {
  const preview = await workbook.render({
    sheetName,
    autoCrop: "all",
    scale: 0.8,
    format: "png",
  });
  await fs.writeFile(
    path.join(previewDir, fileName),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const output = await SpreadsheetFile.exportXlsx(workbook);
const workbookPath = path.join(tableDir, "PATH_CANN_publication_tables_v20.xlsx");
await output.save(workbookPath);
console.log(JSON.stringify({ workbookPath, previewDir }, null, 2));

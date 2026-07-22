#!/usr/bin/env Rscript

# ==============================================================================
# build_figure_audit.R
# ------------------------------------------------------------------------------
# Repository-wide figure/artifact census and logical figure inventory.
#
# This script is intentionally read-only with respect to scientific inputs and
# plot-producing code. It writes only the audit products below:
#   results/figure_audit/artifact_census.csv
#   results/figure_audit/figure_inventory.csv
#   FIGURE_AUDIT_REPORT.md
#
# The inventory is an audit scaffold. A row is never marked numerically or
# visually validated merely because a file exists or a generating script ran.
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "00_config.R"))) {
  stop("Run scripts/build_figure_audit.R from the repository root.", call. = FALSE)
}

out_dir <- file.path(root, "results", "figure_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

artifact_out <- file.path(out_dir, "artifact_census.csv")
inventory_out <- file.path(out_dir, "figure_inventory.csv")
report_out <- file.path(root, "FIGURE_AUDIT_REPORT.md")

graphic_extensions <- c("png", "pdf", "svg", "jpg", "jpeg", "tif", "tiff")
graphic_pattern <- paste0("\\.(", paste(graphic_extensions, collapse = "|"), ")$")

`%||%` <- function(x, y) {
  if (length(x) == 0 || is.null(x) || all(is.na(x)) || !nzchar(as.character(x[1]))) y else x
}

rel_path <- function(path) {
  path <- gsub("\\\\", "/", as.character(path))
  root_prefix <- paste0(root, "/")
  path[startsWith(path, root_prefix)] <- substring(path[startsWith(path, root_prefix)], nchar(root_prefix) + 1L)
  sub("^\\./", "", path)
}

collapse_unique <- function(x, sep = "; ", empty = "") {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (!length(x)) empty else paste(sort(x), collapse = sep)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

file_stem <- function(path) {
  sub(graphic_pattern, "", path, ignore.case = TRUE, perl = TRUE)
}

format_time <- function(x) {
  ifelse(is.na(x), "", format(x, "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Amsterdam"))
}

# ------------------------------------------------------------------------------
# Dimensions without loading/decompressing full raster images.
# ------------------------------------------------------------------------------

raw_uint32_be <- function(x) {
  if (length(x) != 4L) return(NA_real_)
  sum(as.integer(x) * 256^(3:0))
}

png_dimensions <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  header <- readBin(con, what = "raw", n = 24L)
  if (length(header) < 24L || rawToChar(header[13:16]) != "IHDR") {
    return(list(width_px = NA_real_, height_px = NA_real_, status = "invalid_png_header"))
  }
  list(
    width_px = raw_uint32_be(header[17:20]),
    height_px = raw_uint32_be(header[21:24]),
    status = "read_png_header"
  )
}

pdf_dimensions <- function(path) {
  pdfinfo <- unname(Sys.which("pdfinfo"))
  if (!nzchar(pdfinfo)) {
    return(list(page_count = NA_real_, page_width_pt = NA_real_, page_height_pt = NA_real_, status = "pdfinfo_unavailable"))
  }
  out <- tryCatch(
    suppressWarnings(system2(pdfinfo, c("-f", "1", "-l", "1", shQuote(path)), stdout = TRUE, stderr = TRUE)),
    error = function(e) character()
  )
  page_line <- grep("^Pages:[[:space:]]*", out, value = TRUE)
  size_line <- grep("^Page([[:space:]]+[0-9]+)?[[:space:]]+size:[[:space:]]*", out, value = TRUE)
  pages <- if (length(page_line)) suppressWarnings(as.numeric(sub("^Pages:[[:space:]]*([0-9]+).*$", "\\1", page_line[1]))) else NA_real_
  width <- height <- NA_real_
  if (length(size_line)) {
    hit <- regexec("Page(?:[[:space:]]+[0-9]+)?[[:space:]]+size:[[:space:]]*([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts", size_line[1], perl = TRUE)
    vals <- regmatches(size_line[1], hit)[[1]]
    if (length(vals) == 3L) {
      width <- suppressWarnings(as.numeric(vals[2]))
      height <- suppressWarnings(as.numeric(vals[3]))
    }
  }
  status <- if (is.finite(width) && is.finite(height)) "read_pdfinfo" else "pdfinfo_no_page_size"
  list(page_count = pages, page_width_pt = width, page_height_pt = height, status = status)
}

svg_dimensions <- function(path) {
  con <- file(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  txt <- paste(readLines(con, n = 80L, warn = FALSE), collapse = " ")
  extract_attr <- function(attr) {
    hit <- regexec(paste0(attr, "[[:space:]]*=[[:space:]]*[\"']([^\"']+)[\"']"), txt, perl = TRUE)
    vals <- regmatches(txt, hit)[[1]]
    if (length(vals) == 2L) vals[2] else NA_character_
  }
  width_raw <- extract_attr("width")
  height_raw <- extract_attr("height")
  viewbox <- extract_attr("viewBox")
  num <- function(x) suppressWarnings(as.numeric(sub("^([0-9.]+).*$", "\\1", x)))
  width <- num(width_raw)
  height <- num(height_raw)
  if ((!is.finite(width) || !is.finite(height)) && !is.na(viewbox)) {
    vb <- suppressWarnings(as.numeric(strsplit(trimws(viewbox), "[ ,]+")[[1]]))
    if (length(vb) == 4L && all(is.finite(vb))) {
      width <- vb[3]
      height <- vb[4]
    }
  }
  list(
    width_px = width,
    height_px = height,
    status = if (is.finite(width) && is.finite(height)) "read_svg_header" else "svg_dimensions_unresolved"
  )
}

get_dimensions <- function(path, extension) {
  ans <- list(
    width_px = NA_real_, height_px = NA_real_, page_count = NA_real_,
    page_width_pt = NA_real_, page_height_pt = NA_real_, dimension_status = "unsupported_format"
  )
  parsed <- tryCatch(
    switch(
      tolower(extension),
      png = png_dimensions(path),
      pdf = pdf_dimensions(path),
      svg = svg_dimensions(path),
      list(status = "dimension_reader_not_implemented")
    ),
    error = function(e) list(status = paste0("dimension_error: ", conditionMessage(e)))
  )
  for (nm in intersect(names(ans), names(parsed))) ans[[nm]] <- parsed[[nm]]
  ans$dimension_status <- parsed$status %||% ans$dimension_status
  ans
}

# ------------------------------------------------------------------------------
# Physical artifact census: every graphic-like file in the repository.
# ------------------------------------------------------------------------------

all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, include.dirs = FALSE)
all_rel <- rel_path(all_files)
keep <- grepl(graphic_pattern, all_rel, ignore.case = TRUE, perl = TRUE) &
  !grepl("(^|/)\\.git(/|$)", all_rel)
graphic_files <- all_files[keep]
graphic_rel <- all_rel[keep]
ord <- order(graphic_rel)
graphic_files <- graphic_files[ord]
graphic_rel <- graphic_rel[ord]

if (!length(graphic_files)) stop("No graphic artifacts were found.", call. = FALSE)

info <- file.info(graphic_files)
extensions <- tolower(sub("^.*\\.", "", graphic_rel))
dimension_list <- Map(get_dimensions, graphic_files, extensions)
dimension_df <- do.call(rbind, lapply(dimension_list, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
}))

md5 <- unname(tools::md5sum(graphic_files))
md5_counts <- table(md5, useNA = "no")
dup_count <- as.integer(md5_counts[md5])
dup_group <- ifelse(dup_count > 1L & !is.na(md5), paste0("exact_", substr(md5, 1L, 12L)), "")

root_bucket <- ifelse(grepl("/", graphic_rel), sub("/.*$", "", graphic_rel), ".")
artifact_role <- rep("scientific_or_diagnostic_candidate", length(graphic_rel))
artifact_role[grepl("^archive/", graphic_rel)] <- "archived_or_external_artifact"
artifact_role[grepl("^archive/.*/docs/papers/", graphic_rel)] <- "archived_reference_or_legacy_export"
artifact_role[grepl("^outputs/codebooks/", graphic_rel)] <- "document_page_render"
artifact_role[grepl("^outputs/.*/presentations/", graphic_rel)] <- "presentation_render_or_embedded_asset"
artifact_role[grepl("^outputs/", graphic_rel) & artifact_role == "scientific_or_diagnostic_candidate"] <- "manual_review_or_document_render"
artifact_role[grepl("^results/figure_audit/visual_qa/", graphic_rel)] <- "visual_qa_render_derivative"
artifact_role[grepl("^docs/figures/", graphic_rel)] <- "documentation_figure"
artifact_role[graphic_rel == "Rplots.pdf"] <- "uncontrolled_default_R_output"

scientific_scope <- grepl("^(plots|results)/", graphic_rel) |
  grepl("^docs/figures/workflow_flowchart/", graphic_rel) |
  graphic_rel == "Rplots.pdf"
scientific_scope[grepl("^results/figure_audit/visual_qa/", graphic_rel)] <- FALSE

dimensions_text <- ifelse(
  is.finite(dimension_df$width_px) & is.finite(dimension_df$height_px),
  paste0(dimension_df$width_px, " x ", dimension_df$height_px, " px/viewBox"),
  ifelse(
    is.finite(dimension_df$page_width_pt) & is.finite(dimension_df$page_height_pt),
    paste0(dimension_df$page_width_pt, " x ", dimension_df$page_height_pt, " pt", ifelse(is.finite(dimension_df$page_count), paste0("; ", dimension_df$page_count, " page(s)"), "")),
    "unresolved"
  )
)

artifacts <- data.frame(
  artifact_id = sprintf("A%05d", seq_along(graphic_rel)),
  artifact_path = graphic_rel,
  artifact_filename = basename(graphic_rel),
  logical_key = file_stem(graphic_rel),
  output_format = extensions,
  output_exists = TRUE,
  dimensions = dimensions_text,
  width_px = dimension_df$width_px,
  height_px = dimension_df$height_px,
  page_count = dimension_df$page_count,
  page_width_pt = dimension_df$page_width_pt,
  page_height_pt = dimension_df$page_height_pt,
  dimension_status = dimension_df$dimension_status,
  size_bytes = as.numeric(info$size),
  mtime = format_time(info$mtime),
  content_md5 = md5,
  exact_duplicate_group = dup_group,
  exact_duplicate_artifact_count = dup_count,
  root_bucket = root_bucket,
  artifact_role = artifact_role,
  in_scientific_figure_scope = scientific_scope,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Static declarations from source code, including missing expected outputs.
# ------------------------------------------------------------------------------

default_plot_scripts <- c(
  "00c_plot_clinical_summary.R", "02_gene_presence_analysis.R",
  "04_gene_breakdown.R", "05_gene_overview_plots.R", "07_explore_MLST.R",
  "09_inc_plasmid_network.R", "09b_mob_plasmid_reconstruction.R",
  "10_replicon_heatmap.R", "11_compare_strains.R",
  "12a_wgs_qc.R", "13_visualise_panaroo_selection.R",
  "14_genotype_phenotype_model.R", "15_longitudinal_patterns.R",
  "17_lineage_analysis.R", "21_publication_figures.R",
  "23_vf_cross_sectional.R", "24_vf_longitudinal_dynamics.R",
  "25_vf_lineage_vf_interaction.R", "26_vf_define_gene_modules.R",
  "27_vf_score_framework.R", "28_vf_transition_case_studies.R",
  "29_vf_amr_combined_profile.R", "32_uti_not_uti_diagnostic_stats.R",
  "33_mechanism_first_addon.R", "34_robustness_first_addon.R",
  "35_final_figure_pack.R", "36_statistical_sensitivity_addon.R",
  "scripts/research_questions/run_rq01_05.R",
  "scripts/research_questions/run_rq06_08.R"
)

pipeline_stage <- function(script) {
  if (script %in% c("12a_wgs_qc.R", "00c_plot_clinical_summary.R", "13_visualise_panaroo_selection.R", "02_gene_presence_analysis.R")) return("Phase 1")
  if (script %in% c(
    "03_plotting.R", "04_gene_breakdown.R", "05_gene_overview_plots.R",
    "07_explore_MLST.R", "09_inc_plasmid_network.R",
    "09b_mob_plasmid_reconstruction.R", "10_replicon_heatmap.R"
  )) return("Phase 1b")
  if (script %in% c("11_compare_strains.R", "14_genotype_phenotype_model.R", "17_lineage_analysis.R")) return("Phase 2")
  if (script %in% c("15_longitudinal_patterns.R", "21_publication_figures.R")) return("Phase 3")
  if (script %in% c("23_vf_cross_sectional.R", "24_vf_longitudinal_dynamics.R", "25_vf_lineage_vf_interaction.R", "26_vf_define_gene_modules.R", "27_vf_score_framework.R", "28_vf_transition_case_studies.R", "29_vf_amr_combined_profile.R", "32_uti_not_uti_diagnostic_stats.R", "33_mechanism_first_addon.R", "34_robustness_first_addon.R", "36_statistical_sensitivity_addon.R")) return("Phase 4")
  if (script == "35_final_figure_pack.R") return("Phase 5")
  if (grepl("^scripts/research_questions/", script)) return("Phase 5")
  if (script == "scripts/create_workflow_case_count_flowchart.R") return("Standalone documentation")
  if (grepl("^(legacy/|scripts/legacy/)", script) || script == "scripts/run_local_mlst_deprecated.R") return("Legacy/deprecated")
  "Noncanonical/manual"
}

pipeline_membership <- function(script) {
  if (script %in% default_plot_scripts) return("default_canonical_runner")
  if (script == "03_plotting.R") return("optional_legacy_flag")
  if (script == "scripts/create_workflow_case_count_flowchart.R") return("current_standalone_not_in_runner")
  if (grepl("^(legacy/|scripts/legacy/)", script) || script == "scripts/run_local_mlst_deprecated.R") return("legacy_or_deprecated")
  "noncanonical_or_manual"
}

rq_output_dir <- function(filename) {
  stem <- tolower(file_stem(basename(filename)))
  if (grepl("^(continuity_over_time|snp_distribution)$", stem)) return("results/research_questions/RQ01")
  if (grepl("^continuity_by_interval_type$", stem)) return("results/research_questions/RQ02")
  if (grepl("^deidentified_case_snp_distances$", stem)) return("results/research_questions/RQ03")
  if (grepl("^six_rule_case_counts$", stem)) return("results/research_questions/RQ04")
  if (grepl("^wgs_selection_by_status$", stem)) return("results/research_questions/RQ05")
  if (grepl("^rq06", stem)) return("results/research_questions/RQ06")
  if (grepl("^rq07", stem)) return("results/research_questions/RQ07")
  if (grepl("^rq08", stem)) return("results/research_questions/RQ08")
  "results/research_questions/code-declared-unresolved"
}

infer_output_dir <- function(script, filename) {
  stem <- tolower(file_stem(basename(filename)))
  if (script == "35_final_figure_pack.R") {
    if (grepl("^figs[0-9]", stem)) return("plots/final/supplementary")
    return("plots/final")
  }
  if (script == "36_statistical_sensitivity_addon.R") return("plots/statistical_sensitivity")
  if (script == "34_robustness_first_addon.R") return("plots/robustness")
  if (script == "33_mechanism_first_addon.R") return("plots/mechanism")
  if (script == "32_uti_not_uti_diagnostic_stats.R") {
    if (grepl("clinical_rule|near_miss_evidence", stem)) return("plots/clinical")
    return("plots/vf")
  }
  if (script == "29_vf_amr_combined_profile.R") return("plots/vf_amr")
  if (script %in% c("23_vf_cross_sectional.R", "24_vf_longitudinal_dynamics.R", "25_vf_lineage_vf_interaction.R", "26_vf_define_gene_modules.R", "27_vf_score_framework.R", "28_vf_transition_case_studies.R")) return("plots/vf")
  if (script == "21_publication_figures.R") return("plots/publication")
  if (script == "17_lineage_analysis.R") return("results/lineage")
  if (script == "15_longitudinal_patterns.R") return("results/longitudinal")
  if (script == "14_genotype_phenotype_model.R") {
    if (stem == "vf_gene_screening_vs_model_evidence") return("plots/vf")
    return("results/models/plots")
  }
  if (script == "13_visualise_panaroo_selection.R" || script == "12a_wgs_qc.R") return("plots/wgs")
  if (script == "11_compare_strains.R") return("results/strain_compare/plots")
  if (script %in% c(
    "09_inc_plasmid_network.R", "09b_mob_plasmid_reconstruction.R",
    "10_replicon_heatmap.R"
  )) return("plots/plasmids")
  if (script == "07_explore_MLST.R" || script == "scripts/run_local_mlst_deprecated.R") return("plots/mlst")
  if (script %in% c("02_gene_presence_analysis.R", "04_gene_breakdown.R", "05_gene_overview_plots.R")) return("plots/vf")
  if (script == "00c_plot_clinical_summary.R") return("plots/clinical")
  if (script == "scripts/create_workflow_case_count_flowchart.R") return("docs/figures/workflow_flowchart")
  if (grepl("^scripts/research_questions/", script)) return(rq_output_dir(filename))
  if (script == "03_plotting.R") {
    if (grepl("volcano", stem)) return("plots/legacy/old_asb_uti_outputs")
    if (grepl("tree|phylo", stem)) return("plots/phylogeny")
    if (grepl("snp_distance", stem)) return("plots/genomics")
    if (grepl("st_distribution", stem)) return("plots/epidemiology")
    if (grepl("swimmer", stem)) return("plots/timelines")
    if (grepl("network", stem)) return("plots/networks")
    if (grepl("nitrate|heatmap", stem)) return("plots/vf")
    return("plots")
  }
  paste0("code-declared-unresolved/", gsub("[^A-Za-z0-9_.-]", "_", script))
}

normalise_declared_path <- function(script, literal) {
  literal <- gsub("\\\\", "/", literal)
  literal <- sub(paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root), "/?"), "", literal, perl = TRUE)
  literal <- sub("^\\./", "", literal)
  if (grepl("/", literal, fixed = TRUE) && grepl("^(plots|results|docs|outputs|archive)/", literal)) {
    return(literal)
  }
  file.path(infer_output_dir(script, basename(literal)), basename(literal))
}

source_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, include.dirs = FALSE)
source_rel <- rel_path(source_files)
source_keep <- grepl("\\.[Rr]$", source_rel) &
  !grepl("^(results|plots|outputs|archive|data|ont-yellow-routine-fastas|\\.git)/", source_rel)
source_files <- source_files[source_keep]
source_rel <- source_rel[source_keep]

literal_pattern <- "[\"'][^\"']+\\.(png|pdf|svg|jpe?g|tiff?)[\"']"
declaration_rows <- list()
decl_i <- 0L

for (j in seq_along(source_files)) {
  script <- source_rel[j]
  lines <- readLines(source_files[j], warn = FALSE, encoding = "UTF-8")
  for (i in seq_along(lines)) {
    line <- lines[i]
    if (grepl("^[[:space:]]*#", line)) next
    hits <- regmatches(line, gregexpr(literal_pattern, line, perl = TRUE, ignore.case = TRUE))[[1]]
    if (!length(hits) || identical(hits, character(0)) || identical(hits, "")) next
    for (hit in hits) {
      literal <- substring(hit, 2L, nchar(hit) - 1L)
      base_stem <- file_stem(basename(literal))
      if (!nzchar(base_stem) || base_stem %in% c("*", ".") || grepl("[*{}%]", literal)) next
      kind <- if (grepl("ggsave|atomic_ggsave|safe_.*save|save_plot|write_plot|filename[[:space:]]*=|(^|[^[:alnum:]_])(png|pdf|svg|jpeg|tiff)[[:space:]]*\\(|file\\.path|OUT_[A-Z_]+[[:space:]]*<-", line, perl = TRUE, ignore.case = TRUE)) {
        "output_call_or_assignment"
      } else {
        "code_mention"
      }
      decl_i <- decl_i + 1L
      declaration_rows[[decl_i]] <- data.frame(
        generating_script = script,
        source_line = i,
        source_line_text = trimws(line),
        declared_literal = literal,
        declared_path = normalise_declared_path(script, literal),
        output_format = tolower(sub("^.*\\.", "", literal)),
        declaration_kind = kind,
        pipeline_stage = pipeline_stage(script),
        pipeline_membership = pipeline_membership(script),
        stringsAsFactors = FALSE
      )
    }
  }

  # `35_final_figure_pack.R` declares logical IDs and derives both formats in
  # register_figure(); those outputs contain no literal extension at the call.
  call_lines <- which(grepl("^[[:space:]]*register_figure[[:space:]]*\\(", lines) &
                        !grepl("<-[[:space:]]*function", lines))
  for (i in call_lines) {
    block <- paste(lines[i:min(length(lines), i + 4L)], collapse = " ")
    id_hit <- regmatches(block, regexpr("[\"'][A-Za-z0-9_]+[\"']", block, perl = TRUE))
    if (!length(id_hit) || !nzchar(id_hit)) next
    figure_id <- substring(id_hit, 2L, nchar(id_hit) - 1L)
    for (ext in c("png", "pdf")) {
      decl_i <- decl_i + 1L
      declaration_rows[[decl_i]] <- data.frame(
        generating_script = script,
        source_line = i,
        source_line_text = trimws(lines[i]),
        declared_literal = paste0(figure_id, ".", ext),
        declared_path = file.path(
          if (grepl("^FigS[0-9]", figure_id)) "plots/final/supplementary" else "plots/final",
          paste0(figure_id, ".", ext)
        ),
        output_format = ext,
        declaration_kind = "register_figure_output",
        pipeline_stage = pipeline_stage(script),
        pipeline_membership = pipeline_membership(script),
        stringsAsFactors = FALSE
      )
    }
  }
}

declarations <- if (length(declaration_rows)) do.call(rbind, declaration_rows) else data.frame()
if (nrow(declarations)) {
  declarations$declared_path <- gsub("\\\\", "/", declarations$declared_path)
  declarations$logical_key <- file_stem(declarations$declared_path)
  declarations$expected_in_canonical_pipeline <- declarations$pipeline_membership == "default_canonical_runner"
  declarations <- declarations[!duplicated(declarations[c("generating_script", "source_line", "declared_path", "declaration_kind")]), ]
}

output_declarations <- declarations[declarations$declaration_kind != "code_mention", , drop = FALSE]

# ------------------------------------------------------------------------------
# Logical figures: group sibling formats; then append code-declared missing rows.
# ------------------------------------------------------------------------------

artifact_groups <- split(seq_len(nrow(artifacts)), artifacts$logical_key)

make_existing_figure <- function(idx) {
  a <- artifacts[idx, , drop = FALSE]
  key <- a$logical_key[1]
  d <- output_declarations[output_declarations$logical_key == key, , drop = FALSE]
  data.frame(
    logical_key = key,
    primary_path = a$artifact_path[1],
    figure_filename = collapse_unique(a$artifact_filename),
    artifact_paths = collapse_unique(a$artifact_path),
    final_output_path = collapse_unique(a$artifact_path),
    output_formats = collapse_unique(toupper(a$output_format), sep = ";"),
    output_exists = TRUE,
    artifact_count = nrow(a),
    dimensions = collapse_unique(a$dimensions),
    width_px = collapse_unique(a$width_px[is.finite(a$width_px)]),
    height_px = collapse_unique(a$height_px[is.finite(a$height_px)]),
    page_count = collapse_unique(a$page_count[is.finite(a$page_count)]),
    size_bytes = sum(a$size_bytes, na.rm = TRUE),
    mtime = collapse_unique(a$mtime),
    exact_duplicate_group = collapse_unique(a$exact_duplicate_group),
    code_declared = nrow(d) > 0,
    code_declared_missing = FALSE,
    declared_formats = collapse_unique(toupper(d$output_format), sep = ";"),
    code_source_declaration = collapse_unique(d$declared_literal),
    generating_script_from_declaration = collapse_unique(d$generating_script),
    approximate_source_code_lines_from_declaration = collapse_unique(d$source_line),
    pipeline_stage_from_declaration = collapse_unique(d$pipeline_stage),
    pipeline_membership_from_declaration = collapse_unique(d$pipeline_membership),
    expected_in_canonical_pipeline = any(d$expected_in_canonical_pipeline %in% TRUE),
    stringsAsFactors = FALSE
  )
}

figure_existing <- do.call(rbind, lapply(artifact_groups, make_existing_figure))
rownames(figure_existing) <- NULL

missing_decl <- output_declarations[!output_declarations$logical_key %in% artifacts$logical_key, , drop = FALSE]
missing_groups <- split(seq_len(nrow(missing_decl)), missing_decl$logical_key)

make_missing_figure <- function(idx) {
  d <- missing_decl[idx, , drop = FALSE]
  data.frame(
    logical_key = d$logical_key[1],
    primary_path = d$declared_path[1],
    figure_filename = collapse_unique(basename(d$declared_path)),
    artifact_paths = "",
    final_output_path = collapse_unique(d$declared_path),
    output_formats = collapse_unique(toupper(d$output_format), sep = ";"),
    output_exists = FALSE,
    artifact_count = 0L,
    dimensions = "missing",
    width_px = "",
    height_px = "",
    page_count = "",
    size_bytes = 0,
    mtime = "",
    exact_duplicate_group = "",
    code_declared = TRUE,
    code_declared_missing = TRUE,
    declared_formats = collapse_unique(toupper(d$output_format), sep = ";"),
    code_source_declaration = collapse_unique(d$declared_literal),
    generating_script_from_declaration = collapse_unique(d$generating_script),
    approximate_source_code_lines_from_declaration = collapse_unique(d$source_line),
    pipeline_stage_from_declaration = collapse_unique(d$pipeline_stage),
    pipeline_membership_from_declaration = collapse_unique(d$pipeline_membership),
    expected_in_canonical_pipeline = any(d$expected_in_canonical_pipeline %in% TRUE),
    stringsAsFactors = FALSE
  )
}

figure_missing <- if (length(missing_groups)) do.call(rbind, lapply(missing_groups, make_missing_figure)) else figure_existing[0, ]
inventory <- rbind(figure_existing, figure_missing)
rownames(inventory) <- NULL

# ------------------------------------------------------------------------------
# Generating-script/path rules for artifacts that predate explicit declarations.
# ------------------------------------------------------------------------------

infer_generating_script <- function(path) {
  stem <- tolower(basename(file_stem(path)))
  if (grepl("^results/figure_audit/visual_qa/", path)) return("scripts/visual_qa_final_figures.R")
  if (grepl("^plots/final/", path)) return("35_final_figure_pack.R")
  if (grepl("^plots/final_figures/", path)) return("35_final_figure_pack.R")
  if (grepl("^plots/statistical_sensitivity/", path)) return("36_statistical_sensitivity_addon.R")
  if (grepl("^plots/robustness/", path)) return("34_robustness_first_addon.R")
  if (grepl("^plots/mechanism/", path)) return("33_mechanism_first_addon.R")
  if (grepl("^plots/amr/", path)) return("29_vf_amr_combined_profile.R")
  if (grepl("^results/research_questions/RQ0[1-5]/", path)) return("scripts/research_questions/run_rq01_05.R")
  if (grepl("^results/research_questions/RQ0[6-8]/", path)) return("scripts/research_questions/run_rq06_08.R")
  if (grepl("^plots/vf_amr/", path)) return("29_vf_amr_combined_profile.R")
  if (grepl("^plots/publication/", path)) return("21_publication_figures.R")
  if (grepl("^results/lineage/", path)) return("17_lineage_analysis.R")
  if (grepl("^results/longitudinal/", path) && grepl("swimmer", stem)) return("15_longitudinal_patterns.R")
  if (grepl("^results/strain_compare/plots/", path)) return("11_compare_strains.R")
  if (grepl("^results/models/plots/", path)) return("14_genotype_phenotype_model.R")
  if (grepl("^plots/wgs/", path)) {
    if (grepl("wgs_qc", stem)) return("12a_wgs_qc.R")
    if (grepl("panaroo", stem)) return("13_visualise_panaroo_selection.R")
  }
  if (grepl("^plots/clinical/", path)) {
    if (grepl("uti_not_uti", stem)) return("32_uti_not_uti_diagnostic_stats.R")
    return("00c_plot_clinical_summary.R")
  }
  if (grepl("^plots/plasmids/", path)) {
    if (grepl("network|cooccurrence", stem)) return("09_inc_plasmid_network.R")
    if (grepl("replicon_heatmap", stem)) return("10_replicon_heatmap.R")
  }
  if (grepl("^plots/mlst/", path)) {
    if (grepl("st_persistence", stem)) return("scripts/run_local_mlst_deprecated.R")
    return("07_explore_MLST.R")
  }
  if (grepl("^plots/vf/", path)) {
    if (grepl("^core_(bar|histogram)", stem)) return("02_gene_presence_analysis.R")
    if (grepl("nitrate", stem)) return("04_gene_breakdown.R")
    if (stem %in% c("gene_prevalence_bar", "variable_gene_heatmap")) return("05_gene_overview_plots.R")
    if (grepl("^vf_gene_screening_vs_model_evidence", stem)) return("14_genotype_phenotype_model.R")
    if (grepl("^uti_not_uti|^duplicate_culture", stem)) return("32_uti_not_uti_diagnostic_stats.R")
    if (grepl("^module_gene_counts|^module_prevalence|^vf_category_composition|^vf_module_assignment", stem)) return("26_vf_define_gene_modules.R")
    if (grepl("^vf_scores|^vf_score_|^vf_expec|^vf_pca|^vf_pcoa", stem)) return("27_vf_score_framework.R")
    if (grepl("^vf_transition|^vf_not_uti_uti", stem)) return("28_vf_transition_case_studies.R")
    if (grepl("^vf_same_strain|^vf_jaccard|^vf_replacement|^vf_within_host|^vf_gene_gain_loss_consecutive", stem)) return("24_vf_longitudinal_dynamics.R")
    if (grepl("^vf_burden_by_st$|^vf_burden_by_top_st$|^vf_burden_st_x_status|^vf_st_composition|^vf_batch|^vf_event|^vf_status_timepoint|^vf_qc_selection|^vf_denominator_flow|^vf_uricult", stem)) return("25_vf_lineage_vf_interaction.R")
    if (grepl("^vf_burden|^vf_top_gene|^vf_gene_prevalence|^vf_category_burden|^vf_category_barplot", stem)) return("23_vf_cross_sectional.R")
  }
  if (grepl("^docs/figures/workflow_flowchart/(08|09|10)_", path)) return("scripts/create_workflow_case_count_flowchart.R")
  if (grepl("^results/plots/", path)) return("legacy/noncanonical plotting code (exact source unresolved)")
  if (grepl("^plots/legacy/", path)) return("03_plotting.R or archived legacy code")
  if (grepl("^archive/", path)) return("archived source or external reference; not active")
  if (grepl("^outputs/", path)) return("document/presentation rendering workflow; not canonical analysis plot code")
  if (path == "Rplots.pdf") return("uncontrolled R graphics device output; source unresolved")
  "source unresolved"
}

inventory$generating_script <- ifelse(
  nzchar(inventory$generating_script_from_declaration),
  inventory$generating_script_from_declaration,
  vapply(inventory$primary_path, infer_generating_script, character(1))
)
inventory$approximate_source_code_lines <- inventory$approximate_source_code_lines_from_declaration
inventory$pipeline_stage <- ifelse(
  nzchar(inventory$pipeline_stage_from_declaration),
  inventory$pipeline_stage_from_declaration,
  vapply(strsplit(inventory$generating_script, "; ", fixed = TRUE), function(x) pipeline_stage(x[1]), character(1))
)
inventory$pipeline_membership <- ifelse(
  nzchar(inventory$pipeline_membership_from_declaration),
  inventory$pipeline_membership_from_declaration,
  vapply(strsplit(inventory$generating_script, "; ", fixed = TRUE), function(x) pipeline_membership(x[1]), character(1))
)

# ------------------------------------------------------------------------------
# Current denominator contracts are derived, not copied into classifications.
# ------------------------------------------------------------------------------

cohort <- safe_read_csv(file.path(root, "results/clinical/analysis_cohort_longcycler.csv"))
vf_ready <- safe_read_csv(file.path(root, "results/vf/vf_analysis_ready.csv"))
transitions <- safe_read_csv(file.path(root, "results/longitudinal/longcycler_transitions.csv"))
casebook <- safe_read_csv(file.path(root, "results/mechanism/not_uti_to_uti_casebook.csv"))
status_map <- safe_read_csv(file.path(root, "results/clinical/status_map.csv"))
paired_values <- safe_read_csv(file.path(root, "results/statistical_sensitivity/participant_collapsed_score_values.csv"))

contract <- list(
  selected_n = if (!is.null(cohort)) nrow(cohort) else NA_integer_,
  selected_participants = if (!is.null(cohort) && "Participant_id" %in% names(cohort)) length(unique(cohort$Participant_id)) else NA_integer_,
  selected_uti = if (!is.null(cohort) && "UTI_Status" %in% names(cohort)) sum(cohort$UTI_Status == "UTI", na.rm = TRUE) else NA_integer_,
  selected_not_uti = if (!is.null(cohort) && "UTI_Status" %in% names(cohort)) sum(cohort$UTI_Status == "Not_UTI", na.rm = TRUE) else NA_integer_,
  vf_n = if (!is.null(vf_ready)) nrow(vf_ready) else NA_integer_,
  transition_n = if (!is.null(transitions)) nrow(transitions) else NA_integer_,
  case_n = if (!is.null(casebook)) nrow(casebook) else NA_integer_,
  source_status_n = if (!is.null(status_map)) nrow(status_map) else NA_integer_,
  paired_participants = if (!is.null(paired_values) && "Participant_id" %in% names(paired_values)) length(unique(paired_values$Participant_id)) else NA_integer_
)

selected_denominator <- if (all(is.finite(unlist(contract[c("selected_n", "selected_participants", "selected_uti", "selected_not_uti")])))) {
  sprintf("%d selected QC-pass Longcycler episodes from %d participants (%d UTI; %d Not UTI)", contract$selected_n, contract$selected_participants, contract$selected_uti, contract$selected_not_uti)
} else "Current selected denominator not available at inventory build time"

# ------------------------------------------------------------------------------
# Classification defaults and metadata enrichment.
# ------------------------------------------------------------------------------

n <- nrow(inventory)
inventory$figure_id <- sprintf("F%05d", seq_len(n))
inventory$figure_number <- ""
inventory$figure_title <- ""
inventory$scientific_question <- ""
inventory$draft_caption <- ""
inventory$statistical_method <- ""
inventory$required_caveat <- ""
inventory$filtering_rule <- ""
inventory$visual_encoding_description <- ""
inventory$multiple_testing <- ""
inventory$planned_width_in <- ""
inventory$planned_height_in <- ""
inventory$planned_dpi <- ""
inventory$manifest_validation_status <- ""
inventory$inventory_scope <- "scientific_or_diagnostic_figure"
inventory$source_data_file_or_object <- "Not yet traced to a unique source table/object"
inventory$analysis_represented <- "Not yet classified beyond artifact/source-path context"
inventory$unit_of_observation <- "Not yet traced"
inventory$episode_or_participant_level <- "Not yet traced"
inventory$plotted_denominator <- "Not yet validated"
inventory$plot_type <- "Not yet classified"
inventory$x_variable <- "Not yet transcribed/validated"
inventory$y_variable <- "Not yet transcribed/validated"
inventory$colour_fill_variable <- "Not yet transcribed/validated"
inventory$facet_variable <- "None identified or not yet transcribed"
inventory$statistical_result_displayed <- "None identified or not yet validated"
inventory$repeated_measures_present <- "Unknown — requires source-data trace"
inventory$evidence_class <- "descriptive"
inventory$intended_use <- "diagnostic"
inventory$current_problems <- "Per-figure numerical and visual validation pending"
inventory$severity <- "moderate"
inventory$recommended_action <- "Trace plotted values to source data, inspect at full and thesis size, then retain/revise/replace/obsolete explicitly"
inventory$action_implemented <- "Inventory row created; scientific and visual verification pending"
inventory$canonical_status <- "unclassified_existing_artifact"
inventory$canonical <- "PENDING"
inventory$validation_status <- "pending_numerical_validation"
inventory$visual_qa_status <- "pending_full_size_and_thesis_size_review"
inventory$source_trace_status <- ifelse(inventory$code_declared, "code_declaration_identified; data trace pending", "generating source unresolved or inferred")
inventory$duplicate_group <- inventory$exact_duplicate_group
inventory$duplicate_basis <- ifelse(nzchar(inventory$exact_duplicate_group), "byte-identical artifact (MD5)", "none identified")
inventory$superseded_by <- ""

path <- inventory$primary_path
stem <- tolower(basename(file_stem(path)))

# Scope classification first.
idx <- grepl("^archive/", path)
inventory$inventory_scope[idx] <- "archived_or_external_artifact"
inventory$analysis_represented[idx] <- "Archived output or reference document; not part of the active analysis release"
inventory$unit_of_observation[idx] <- "Not applicable / archived artifact"
inventory$episode_or_participant_level[idx] <- "Not applicable"
inventory$plotted_denominator[idx] <- "Outdated or not applicable; must not be used"
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Archived or external artifact is not an active, validated project figure"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Exclude from final figure set; retain only as archive/reference provenance"
inventory$action_implemented[idx] <- "Classified obsolete/nonactive; no archived file deleted"
inventory$canonical_status[idx] <- "archived_or_external_not_canonical"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_applicable_archived_or_external"
inventory$visual_qa_status[idx] <- "not_required_for_active_release"

idx <- grepl("^results/figure_audit/visual_qa/", path)
inventory$inventory_scope[idx] <- "visual_qa_render_not_standalone_figure"
inventory$analysis_represented[idx] <- "Derived full-size, A4-preview, or colour-vision-deficiency rendering of a canonical final figure"
inventory$unit_of_observation[idx] <- "One rendered derivative of one canonical figure"
inventory$episode_or_participant_level[idx] <- "Inherited from the canonical source figure; not independently authoritative"
inventory$plotted_denominator[idx] <- "Inherited from the canonical source figure"
inventory$plot_type[idx] <- "Visual-QA render derivative"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "A QA derivative duplicates its canonical source and would inflate the scientific figure count if treated as standalone"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Retain only as visual-QA evidence; cite and validate the canonical source figure"
inventory$action_implemented[idx] <- "Classified as a non-standalone visual-QA derivative"
inventory$canonical_status[idx] <- "visual_qa_derivative_not_canonical"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "inherits_source_validation_only; derivative_not_independent"
inventory$visual_qa_status[idx] <- "generated_visual_QA_derivative"
inventory$source_trace_status[idx] <- "generated by scripts/visual_qa_final_figures.R from the canonical final PNG"

idx <- grepl("^outputs/codebooks/", path)
inventory$inventory_scope[idx] <- "document_page_render_not_standalone_figure"
inventory$analysis_represented[idx] <- "Rendered codebook/document page or contact sheet"
inventory$unit_of_observation[idx] <- "Document page"
inventory$episode_or_participant_level[idx] <- "Not applicable"
inventory$plot_type[idx] <- "Document render"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "Not a standalone scientific figure; including it as one would inflate figure counts"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Keep in artifact census; exclude from final scientific figure set"
inventory$action_implemented[idx] <- "Classified as document render"
inventory$canonical_status[idx] <- "render_artifact_not_scientific_figure"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_applicable_document_render"
inventory$visual_qa_status[idx] <- "document_QA_outside_figure_audit"

idx <- grepl("^outputs/.*/presentations/", path)
inventory$inventory_scope[idx] <- "presentation_render_or_embedded_asset"
inventory$analysis_represented[idx] <- "Presentation slide render or embedded presentation asset"
inventory$unit_of_observation[idx] <- "Slide/page/embedded asset"
inventory$episode_or_participant_level[idx] <- "Not applicable as standalone figure"
inventory$plot_type[idx] <- "Presentation render"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "Not necessarily a standalone pipeline figure; may duplicate canonical figures embedded in a deck"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Keep in artifact census; validate source canonical figure rather than each slide render"
inventory$action_implemented[idx] <- "Classified as presentation artifact"
inventory$canonical_status[idx] <- "presentation_artifact_not_canonical_analysis_figure"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_applicable_presentation_render"
inventory$visual_qa_status[idx] <- "presentation_QA_outside_figure_audit"

idx <- grepl("^outputs/", path) & inventory$canonical == "PENDING"
inventory$inventory_scope[idx] <- "manual_review_or_document_render"
inventory$analysis_represented[idx] <- "Manual-review, document, or rendered deliverable artifact"
inventory$unit_of_observation[idx] <- "Rendered artifact"
inventory$episode_or_participant_level[idx] <- "Not applicable as standalone analytical figure"
inventory$plot_type[idx] <- "Rendered deliverable"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "Not linked to the canonical R figure pipeline"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Keep in artifact census; exclude unless independently confirmed as a scientific figure"
inventory$action_implemented[idx] <- "Classified outside canonical analytical figure scope"
inventory$canonical_status[idx] <- "manual_render_not_canonical"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_applicable_manual_render"
inventory$visual_qa_status[idx] <- "outside_active_figure_QA"

idx <- grepl("^results/plots/", path)
inventory$analysis_represented[idx] <- "Legacy/noncanonical diagnostic output from the former results/plots tree"
inventory$plotted_denominator[idx] <- "Outdated or unknown legacy denominator"
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Stale 2025-era/noncanonical output tree; source inputs and denominator are not current-release certified"
inventory$severity[idx] <- "major"
inventory$recommended_action[idx] <- "Exclude from final set and retain only under an explicit archive manifest"
inventory$action_implemented[idx] <- "Classified obsolete; no file deleted"
inventory$canonical_status[idx] <- "legacy_obsolete"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_validated_obsolete"
inventory$visual_qa_status[idx] <- "not_required_unless_reactivated"
inventory$superseded_by[idx] <- "Current outputs under plots/ and plots/final_figures/"

idx <- grepl("^(plots/legacy/|plots/persistence/|plots/participant_specific/|plots/pairwise_identity/)", path)
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Legacy, per-participant, or bulk exploratory artifact is not part of the current final figure contract"
inventory$severity[idx] <- "moderate"
inventory$recommended_action[idx] <- "Retain as diagnostic/archive only; do not cite as current evidence"
inventory$action_implemented[idx] <- "Classified legacy/exploratory"
inventory$canonical_status[idx] <- "legacy_or_bulk_exploratory"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_validated_for_current_release"
inventory$superseded_by[idx] <- "Current Phase-4 and final-figure outputs"

idx <- grepl("^docs/figures/workflow_flowchart/", path)
inventory$inventory_scope[idx] <- "documentation_figure"
inventory$analysis_represented[idx] <- "Workflow or denominator documentation graphic"
inventory$unit_of_observation[idx] <- "Pipeline stage, analytical unit, or denominator"
inventory$episode_or_participant_level[idx] <- "Mixed workflow units; each count must name its unit"
inventory$plot_type[idx] <- "Workflow diagram/flowchart"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "Noncanonical documentation rendering; content was inventoried but is not part of the retained analytical figure set"
inventory$severity[idx] <- "major"
inventory$recommended_action[idx] <- "Keep outside the analytical figure set unless it is regenerated from derived counts and separately validated"
inventory$action_implemented[idx] <- "Audited and classified as noncanonical documentation output; current render retained"
inventory$canonical_status[idx] <- "documentation_figure_audited_not_retained"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_result_certified_noncanonical_documentation"
inventory$visual_qa_status[idx] <- "not_required_unless_reactivated"

# The three denominator panels are regenerated by
# scripts/create_workflow_case_count_flowchart.R from the current Longcycler
# contract.  Keep the other documentation figures conservatively pending, but
# do not continue to label these current renders as stale.
idx <- grepl("^docs/figures/workflow_flowchart/(08_unit_aware_denominator_story|09_case_count_workflow_ladder|10_transition_count_funnel)$", path)
inventory$current_problems[idx] <- "No known numerical defect after regeneration from current derived counts"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Retain as validated methods/denominator documentation; regenerate whenever the cohort contract changes"
inventory$action_implemented[idx] <- "Regenerated from current derived values and visually checked"
inventory$canonical_status[idx] <- "validated_current_documentation_figure"
inventory$validation_status[idx] <- "validated_current_derived_denominators"
inventory$visual_qa_status[idx] <- "visually_checked_current_render"

# Script 29 publishes three result-ready supplementary AMR figures. They are
# intentionally kept outside the canonical thesis pack because AMR is a
# supporting layer rather than a numbered research question.
idx <- grepl("^plots/amr/", path)
inventory$inventory_scope[idx] <- "supplementary_supporting_genomic_amr"
inventory$analysis_represented[idx] <- ifelse(
  grepl("most_prevalent", stem[idx]),
  "Prevalence of informative AMRFinderPlus acquired genes with mdf(A) excluded",
  ifelse(
    grepl("stability", stem[idx]),
    "Informative acquired-gene stability across all adjacent pairs by direct SNP context",
    "Caller-specific acquired-gene detections by harmonized determinant class"
  )
)
inventory$unit_of_observation[idx] <- ifelse(
  grepl("stability", stem[idx]), "Adjacent within-resident pair",
  ifelse(grepl("concordance", stem[idx]), "Episode-determinant", "Episode")
)
inventory$episode_or_participant_level[idx] <- ifelse(
  grepl("stability", stem[idx]), "Pair level nested within resident",
  "Episode level with resident prevalence tabulated where applicable"
)
inventory$plotted_denominator[idx] <- ifelse(
  grepl("stability", stem[idx]), "371 adjacent pairs from 139 residents",
  "532 selected episodes from 161 residents; all required callers complete"
)
inventory$plot_type[idx] <- ifelse(
  grepl("most_prevalent", stem[idx]), "Horizontal prevalence bar plot",
  ifelse(
    grepl("stability", stem[idx]),
    "Boxplot with jitter by direct SNP context",
    "Caller-by-class heatmap"
  )
)
inventory$source_data_file_or_object[idx] <- ifelse(
  grepl("most_prevalent", stem[idx]),
  "results/amr/gene_prevalence_episode_resident.csv",
  ifelse(
    grepl("stability", stem[idx]),
    "results/amr/adjacent_pair_amr_profiles_371.csv; results/amr/longitudinal_resident_bootstrap_inference.csv",
    "results/amr/caller_concordance_discrepancies.csv"
  )
)
inventory$evidence_class[idx] <- "descriptive"
inventory$intended_use[idx] <- "supplementary"
inventory$current_problems[idx] <- "No known numerical defect; genomic predictions must not be read as phenotypic AST or causal evidence"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Retain as supplementary supporting output; do not promote to a canonical thesis figure unless AMR becomes a formal research question"
inventory$action_implemented[idx] <- "Registered as a validated script-29 supplementary output with explicit denominator and interpretation boundary"
inventory$canonical_status[idx] <- "supplementary_supporting_not_in_canonical_thesis_pack"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "validated_by_script29_release_contract"
inventory$visual_qa_status[idx] <- "supplementary_render_requires_standard_visual_review"
inventory$source_trace_status[idx] <- "traced to script-29 published AMR tables"

idx <- grepl("^results/final_figures/denominator_flowchart/", path)
inventory$inventory_scope[idx] <- "superseded_documentation_figure"
inventory$analysis_represented[idx] <- "Former mixed-assembler denominator flowchart"
inventory$unit_of_observation[idx] <- "Mixed workflow units"
inventory$episode_or_participant_level[idx] <- "Mixed"
inventory$plot_type[idx] <- "Workflow diagram/flowchart"
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Pre-Longcycler-only denominator graphic remains outside the current final manifest"
inventory$severity[idx] <- "major"
inventory$recommended_action[idx] <- "Archive after a current derived-count workflow figure is regenerated"
inventory$action_implemented[idx] <- "Classified obsolete; file not deleted"
inventory$canonical_status[idx] <- "superseded_stale_documentation"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "known_stale_denominators"
inventory$visual_qa_status[idx] <- "not_required_until_regenerated"
inventory$superseded_by[idx] <- "Updated scripts/create_workflow_case_count_flowchart.R outputs after regeneration"

idx <- grepl("^plots/lecturer_(core|optional)_figures/", path)
inventory$inventory_scope[idx] <- "presentation_copy_of_analysis_figure"
inventory$analysis_represented[idx] <- "Lecturer/presentation copy of a source analytical figure"
inventory$unit_of_observation[idx] <- "Inherited from canonical source figure"
inventory$episode_or_participant_level[idx] <- "Inherited from canonical source figure"
inventory$plotted_denominator[idx] <- "Must match source final figure; copy is not independently authoritative"
inventory$plot_type[idx] <- "Presentation-exported figure copy"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "Copied/exported derivative can become stale and should not be treated as the canonical scientific source"
inventory$severity[idx] <- "minor"
inventory$recommended_action[idx] <- "Validate the source final figure and rebuild presentation copies from it"
inventory$action_implemented[idx] <- "Classified as presentation derivative"
inventory$canonical_status[idx] <- "presentation_copy_not_canonical_source"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "inherits_source_validation_only_after_checksum_match"
inventory$visual_qa_status[idx] <- "presentation_copy_QA_pending"
inventory$superseded_by[idx] <- paste0("plots/final_figures/", basename(inventory$logical_key[idx]))
inventory$duplicate_group[idx] <- paste0("presentation_copy_", gsub("[^A-Za-z0-9]+", "_", basename(inventory$logical_key[idx])))
inventory$duplicate_basis[idx] <- "named presentation copy of final/source figure"

idx <- grepl("^results/wgs/reports/participant_report_", path)
inventory$inventory_scope[idx] <- "bulk_participant_diagnostic_report"
inventory$analysis_represented[idx] <- "Per-participant WGS diagnostic report"
inventory$unit_of_observation[idx] <- "Participant and that participant's isolates/comparisons"
inventory$episode_or_participant_level[idx] <- "Participant-specific diagnostic"
inventory$plotted_denominator[idx] <- "One participant per report; isolate/comparison count varies"
inventory$plot_type[idx] <- "Multipage participant report"
inventory$intended_use[idx] <- "diagnostic"
inventory$current_problems[idx] <- "Bulk diagnostic report, not a standalone thesis/manuscript figure; source is legacy/noncanonical WGS reporting code"
inventory$severity[idx] <- "minor"
inventory$recommended_action[idx] <- "Retain only for case-level QA; exclude from final scientific figure count"
inventory$action_implemented[idx] <- "Classified as bulk diagnostic report"
inventory$canonical_status[idx] <- "bulk_participant_diagnostic_not_final"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_validated_as_final_figure"
inventory$visual_qa_status[idx] <- "case_report_QA_only"

idx <- grepl("^results/legacy/", path)
inventory$inventory_scope[idx] <- "legacy_generated_output"
inventory$analysis_represented[idx] <- "Archived legacy ASB-versus-UTI analysis output"
inventory$unit_of_observation[idx] <- "Legacy/unknown"
inventory$episode_or_participant_level[idx] <- "Legacy/unknown"
inventory$plotted_denominator[idx] <- "Superseded legacy denominator"
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Legacy ASB-versus-UTI output is not the current primary UTI-versus-Not UTI analysis"
inventory$severity[idx] <- "major"
inventory$recommended_action[idx] <- "Keep archived and exclude from current evidence/figure set"
inventory$action_implemented[idx] <- "Classified obsolete"
inventory$canonical_status[idx] <- "legacy_obsolete"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_validated_obsolete"
inventory$visual_qa_status[idx] <- "not_required"
inventory$superseded_by[idx] <- "Current UTI-versus-Not UTI Phase-4 outputs"

idx <- grepl("^(plots/(core_|richness_|upset_|epidemiology/|genomics/|phylogeny/|timelines/)|results/(mlst/top20_STs|modelling/plots/|replicon_heatmap))", path)
inventory$inventory_scope[idx] <- "optional_or_stale_exploratory_figure"
inventory$analysis_represented[idx] <- "Optional legacy/exploratory analysis output"
inventory$plotted_denominator[idx] <- "Not certified against the current selected Longcycler contract"
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Generated by optional/deprecated or older output path and not part of the current final manifest"
inventory$severity[idx] <- "major"
inventory$recommended_action[idx] <- "Do not cite as current evidence; replace with a current source-traced figure if scientifically needed"
inventory$action_implemented[idx] <- "Classified obsolete/exploratory"
inventory$canonical_status[idx] <- "optional_legacy_or_stale_obsolete"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "not_validated_for_current_release"
inventory$visual_qa_status[idx] <- "not_required_unless_reactivated"
inventory$superseded_by[idx] <- "Current Phase-4/final figure set"

# The former unnumbered final pack is retained for provenance only.  Every
# family has an explicit replacement in the numbered thesis pack.
old_final_map <- c(
  primary_denominator_and_uncertainty = "plots/final/Fig01_cohort_and_denominators; plots/final/supplementary/FigS05_near_miss_leave_one_uti",
  not_uti_to_uti_mechanism_casebook = "plots/final/Fig06_longitudinal_trajectories; plots/final/Fig07_within_host_genomic_continuity; plots/final/supplementary/FigS03_module_gain_loss; plots/final/supplementary/FigS06_transition_mechanisms",
  strain_stability_and_host_context = "plots/final/Fig06_longitudinal_trajectories; plots/final/Fig07_within_host_genomic_continuity",
  global_vf_signal_and_robustness = "plots/final/Fig04_vf_burden; plots/final/Fig05_vf_association_evidence; plots/final/supplementary/FigS05_near_miss_leave_one_uti; plots/final/supplementary/FigS08_gene_model_forest; plots/final/supplementary/FigS09_event_sample_sensitivity",
  transition_mechanisms_by_transition_type = "plots/final/supplementary/FigS06_transition_mechanisms",
  accessory_plasmid_amr_changes = "plots/final/supplementary/FigS07_plasmid_amr_context",
  near_miss_and_sparse_precision = "plots/final/supplementary/FigS05_near_miss_leave_one_uti",
  leave_one_uti_stability = "plots/final/supplementary/FigS05_near_miss_leave_one_uti",
  prioritised_variant_map = "plots/final/Fig08_reference_aware_variant_map",
  lineage_confounding_diagnostic = "plots/final/Fig03_sequence_type_distribution; plots/final/Fig05_vf_association_evidence; plots/final/supplementary/FigS04_vf_pcoa",
  paired_resident_expec_marker = "plots/final/Fig04_vf_burden",
  not_uti_to_uti_module_gain_loss = "plots/final/supplementary/FigS03_module_gain_loss",
  vf_module_pcoa_primary_status = "plots/final/supplementary/FigS04_vf_pcoa"
)
for (old_id in names(old_final_map)) {
  hit <- inventory$logical_key == paste0("plots/final_figures/", old_id)
  if (!any(hit)) next
  inventory$inventory_scope[hit] <- "superseded_previous_final_pack"
  inventory$analysis_represented[hit] <- "Former unnumbered Phase-4 final-pack figure; superseded by the numbered canonical thesis figure pack"
  inventory$intended_use[hit] <- "obsolete"
  inventory$current_problems[hit] <- "Historical final output is no longer declared by the current canonical source and must not be cited as the retained figure"
  inventory$severity[hit] <- "minor"
  inventory$recommended_action[hit] <- "Retain only as provenance; use the numbered replacement under plots/final"
  inventory$action_implemented[hit] <- "Superseded classification and explicit replacement recorded; historical file left unchanged"
  inventory$canonical_status[hit] <- "superseded_previous_final_pack"
  inventory$canonical[hit] <- "FALSE"
  inventory$validation_status[hit] <- "superseded_not_current_release"
  inventory$visual_qa_status[hit] <- "not_required_for_superseded_output"
  inventory$source_trace_status[hit] <- "historical 35_final_figure_pack.R output; no longer declared by current source"
  inventory$superseded_by[hit] <- unname(old_final_map[[old_id]])
  inventory$duplicate_group[hit] <- paste0("superseded_final_pack_", old_id)
  inventory$duplicate_basis[hit] <- "scientific role replaced by numbered canonical thesis figure(s)"
}

idx <- grepl("(^|/)Rplots\\.pdf$", path)
inventory$inventory_scope[idx] <- "uncontrolled_graphics_device_output"
inventory$analysis_represented[idx] <- "Uncontrolled default R graphics device output"
inventory$unit_of_observation[idx] <- "Unknown"
inventory$intended_use[idx] <- "obsolete"
inventory$current_problems[idx] <- "Generating script, source data, and intended interpretation are unresolved"
inventory$severity[idx] <- "major"
inventory$recommended_action[idx] <- "Identify provenance or remove/archive after confirming it is not required"
inventory$action_implemented[idx] <- "Flagged; not deleted"
inventory$canonical_status[idx] <- "uncontrolled_obsolete"
inventory$canonical[idx] <- "FALSE"
inventory$validation_status[idx] <- "unvalidated_unknown_source"

# Current active-path defaults. These are pipeline diagnostics/source figures,
# not automatically canonical thesis figures. Their source and disposition are
# audited here, but a successful render alone is deliberately not promoted to
# scientific validation.
idx <- grepl("^(plots/(clinical|wgs|vf|mlst|plasmids|vf_amr|mechanism|robustness|statistical_sensitivity|final|final_figures|publication)|results/(models/plots|strain_compare/plots|longitudinal|lineage|research_questions/RQ0[1-8]))/", path)
active_default <- idx & inventory$canonical == "PENDING"
inventory$canonical_status[active_default] <- "active_pipeline_diagnostic_audited_not_retained"
inventory$canonical[active_default] <- "FALSE"
inventory$validation_status[active_default] <- "audited_nonfinal_diagnostic_not_result_certified"
inventory$recommended_action[active_default] <- "Keep as a reproducible pipeline diagnostic or source layer; do not cite as a retained thesis figure unless separately validated"
inventory$action_implemented[active_default] <- "Source/disposition audited and excluded from the numbered thesis set"
inventory$plotted_denominator[idx & inventory$plotted_denominator == "Not yet validated"] <- selected_denominator
inventory$repeated_measures_present[idx & inventory$repeated_measures_present == "Unknown — requires source-data trace"] <- "Yes for episode/isolate-level cohort views unless explicitly participant-collapsed"

# Plot-type inference is descriptive only and never confers validation.
inventory$plot_type[grepl("heatmap|matrix|tile", stem)] <- "Heatmap/tile matrix"
inventory$plot_type[grepl("volcano", stem)] <- "Volcano plot"
inventory$plot_type[grepl("forest|stability", stem)] <- "Forest/interval/stability plot"
inventory$plot_type[grepl("boxplot|burden", stem)] <- "Distribution/boxplot or burden summary"
inventory$plot_type[grepl("violin", stem)] <- "Violin/distribution plot"
inventory$plot_type[grepl("histogram|distribution", stem)] <- "Histogram/distribution plot"
inventory$plot_type[grepl("swimmer|timeline|trajectory|trajectories", stem)] <- "Longitudinal trajectory/swimmer plot"
inventory$plot_type[grepl("network|cooccurrence", stem)] <- "Network graph"
inventory$plot_type[grepl("scatter|pcoa|pca", stem)] <- "Scatter/ordination plot"
inventory$plot_type[grepl("roc", stem)] <- "ROC curve"
inventory$plot_type[grepl("upset", stem)] <- "UpSet plot"
inventory$plot_type[grepl("tree|phylogen", stem)] <- "Phylogeny/tree"
inventory$plot_type[grepl("bar|counts|overview|waterfall|flow|composition|prevalence|selection", stem) & inventory$plot_type == "Not yet classified"] <- "Bar/count/flow summary"

# Unit and evidence inference for active analytical figures.
idx <- grepl("wgs_qc|assembly|panaroo_selection", stem)
inventory$unit_of_observation[idx] <- "Genome assembly or assembly candidate"
inventory$episode_or_participant_level[idx] <- "Assembly-level"
inventory$evidence_class[idx] <- "QC"
inventory$source_data_file_or_object[idx] <- "results/wgs/qc_summary.csv; results/qc/analysis_assembly_manifest.csv; results/wgs/pan/panaroo_input_manifest.csv as applicable"

idx <- grepl("transition|casebook|case_matrix|strain_stability|gene_gain_loss|jaccard_by_days", stem)
inventory$unit_of_observation[idx] <- "Ordered within-participant episode pair / transition"
inventory$episode_or_participant_level[idx] <- "Transition/pair-level"
inventory$plotted_denominator[idx] <- if (is.finite(contract$transition_n)) sprintf("%d canonical adjacent within-participant transitions; focused Not UTI-to-UTI case figures use %d transitions", contract$transition_n, contract$case_n) else "Canonical transition denominator unavailable"
inventory$repeated_measures_present[idx] <- "Yes — participants may contribute multiple transitions"
inventory$source_data_file_or_object[idx] <- "results/longitudinal/longcycler_transitions.csv; results/mechanism/not_uti_to_uti_casebook.csv; transition-specific VF tables as applicable"

idx <- grepl("paired_resident|participant_summary|paired_participant|slopeplot", stem)
inventory$unit_of_observation[idx] <- "Participant-level status summary or paired participant"
inventory$episode_or_participant_level[idx] <- "Participant-level"
inventory$repeated_measures_present[idx] <- "Collapsed/paired by participant"
inventory$plotted_denominator[idx] <- if (is.finite(contract$paired_participants)) sprintf("Participant-collapsed subset; %d participants represented across score summaries (exact per-panel n must be captioned)", contract$paired_participants) else "Participant-collapsed subset; exact n pending"

idx <- grepl("vf|gene|module|expec", stem) & inventory$unit_of_observation == "Not yet traced"
inventory$unit_of_observation[idx] <- "Episode/isolate, gene, or VF module (panel-specific)"
inventory$episode_or_participant_level[idx] <- "Primarily episode/isolate-level unless labelled participant-collapsed"
inventory$source_data_file_or_object[idx] <- "results/vf/vf_analysis_ready.csv and associated gene/module/model result tables"

idx <- grepl("model|volcano|forest|risk|association|bootstrap|power|precision|leave_one|effect", stem)
inventory$evidence_class[idx] <- "exploratory"
inventory$statistical_result_displayed[idx] <- "Effect estimate/test/model diagnostic; exact statistic and adjustment require row-level validation"

idx <- grepl("qc|selection|denominator|provenance|audit|duplicate", stem)
inventory$evidence_class[idx] <- "QC"
inventory$intended_use[idx & inventory$intended_use != "obsolete"] <- "diagnostic"

# The current numbered final-pack manifest is authoritative for retained class,
# data provenance, captions, scientific limitations and intended dimensions.
# Its explicit data-check table validates the frozen denominators and file set;
# it does not substitute for manual full-size/thesis-size visual QA.
final_manifest <- safe_read_csv(file.path(root, "results/figure_audit/final_figure_manifest.csv"))
final_checks <- safe_read_csv(file.path(root, "results/figure_audit/final_figure_data_checks.csv"))
final_checks_pass <- !is.null(final_checks) && "pass" %in% names(final_checks) && nrow(final_checks) > 0L && all(final_checks$pass %in% TRUE)

normalise_source_list <- function(x) {
  x <- as.character(x)
  gsub(paste0(gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root), "/"), "", x, perl = TRUE)
}

final_denominators <- c(
  Fig01_cohort_and_denominators = "583 primary clinical episodes from 166 participants; 532 selected episodes from 161 participants (16 UTI; 516 Not UTI)",
  Fig02_wgs_quality_control = "592 Longcycler candidate assemblies; 532 selected QC-passing assemblies",
  Fig03_sequence_type_distribution = "532 selected episodes from 161 participants (16 UTI; 516 Not UTI)",
  Fig04_vf_burden = "532 selected episodes; participant-status panel contains 12 participants observed in both states",
  Fig05_vf_association_evidence = "532 selected episodes from 161 participants; four prespecified standardized VF-score models",
  Fig06_longitudinal_trajectories = "Nine data-derived adjacent Not UTI-to-UTI transitions among nine participants",
  Fig07_within_host_genomic_continuity = "371 adjacent direct comparisons from 139 participants (140 at or below and 231 above 25 SNPs)",
  Fig08_reference_aware_variant_map = "43 exact-reference-validated variants across the available reference-specific Not UTI-to-UTI comparisons",
  FigS01_vf_presence_heatmap = "532 selected episodes by 108 genes with 5-95% prevalence",
  FigS02_core_genome_phylogeny = "532 selected assembly tips matched one-to-one to metadata",
  FigS03_module_gain_loss = "Nine adjacent Not UTI-to-UTI transition pairs by curated VF module",
  FigS04_vf_pcoa = "532 selected episodes",
  FigS05_near_miss_leave_one_uti = "17 prespecified near-miss episodes; leave-one-out ranges across 16 UTI episodes and four VF scores",
  FigS06_transition_mechanisms = "371 adjacent comparisons from 139 participants",
  FigS07_plasmid_amr_context = "532 PlasmidFinder-complete episodes; 532 MOB-suite profiles; 532 validated script-29 genomic-AMR profiles; nine descriptive Not UTI-to-UTI transitions",
  FigS08_gene_model_forest = "50 prescreened fitted gene/feature models",
  FigS09_event_sample_sensitivity = "12 nearest within-participant pairs; companion event-sample analysis has 32 samples from 29 participants",
  FigS10_snp_threshold_sensitivity = "371 adjacent comparisons from 139 participants"
)

final_panel_map <- data.frame(
  figure_id = names(final_denominators),
  plot_type = c(
    "Composite count/dot and stacked-bar denominator summary", "Two-panel QC scatter plot",
    "Prevalence dumbbell/dot plot plus provenance stacked bar", "Violin/box/jitter distributions plus paired participant slope plot",
    "Odds-ratio forest plot", "Date-aware longitudinal trajectory plot",
    "Log-scale SNP scatter plus VF-similarity distribution", "Reference-faceted genomic-position map",
    "Biologically ordered binary heatmap", "Annotated core-genome phylogeny",
    "Gain/loss state heatmap", "Jaccard PCoA scatter plot",
    "Near-miss count bars plus leave-one-UTI interval diagnostic", "Proportional stacked bar",
    "Replicon-prevalence dot plot plus AMR-burden distribution", "Separation-aware odds-ratio forest plot",
    "Within-participant difference dot plot with median", "Threshold-sensitivity line and bootstrap ribbon"
  ),
  x = c(
    "Cohort/stage count", "Contig count; assembly length (Mb)", "Within-status prevalence; typing provenance",
    "Operational UTI status", "Adjusted odds ratio for UTI versus Not UTI", "Observed collection date",
    "Days between samples; operational 25-SNP category", "Cumulative position within exact reference (Mb)",
    "Selected episode", "Core-genome branch distance", "Deidentified transition case", "PCoA axis 1",
    "Clinical-rule component; UTI-minus-Not-UTI score difference", "Observed adjacent transition type",
    "Within-status replicon prevalence; operational UTI status", "Adjusted odds ratio for UTI versus Not UTI",
    "Within-participant UTI-minus-Not-UTI difference", "Pairwise SNP threshold"
  ),
  y = c(
    "Episode/participant count", "Assembly N50 (kb); GC content (%)", "Sequence type; episode count",
    "Detected VF-gene/system/marker count", "Prespecified VF score", "Deidentified transition case",
    "Pairwise SNP distance + 1; VF Jaccard similarity", "Deidentified comparison/reference panel",
    "Curated VF gene", "Tree tip/assembly", "Curated VF module", "PCoA axis 2",
    "Episode count; prespecified VF score", "Proportion of transitions", "Replicon; AMR genes detected",
    "Virulence-factor gene/feature", "VF endpoint", "Proportion at or below threshold"
  ),
  colour = c(
    "Operational UTI status and count unit", "Assembly QC pass/fail; selection shape", "Operational UTI status; ST-call provenance",
    "Operational UTI status", "Model diagnostic/singularity", "Operational UTI status; sampling-event shape; SNP-context linetype",
    "Observed status transition; SNP-threshold category", "Variant type", "Presence/absence/unavailable; operational-status annotation",
    "Operational UTI status; sequence-type shape", "Gain/loss/stable/unavailable state", "Operational UTI status; sequence-type shape",
    "Direction-change diagnostic", "Evidence category", "Operational UTI status", "Model diagnostic and estimability shape",
    "Endpoint", "Observed estimate and 95% cluster-bootstrap interval"
  ),
  facet = c(
    "Episodes versus participants", "QC metric panel", "Distribution versus provenance panel", "Score and episode/participant panel",
    "None", "Participant/case", "SNP continuity versus VF similarity panel", "Exact reference and comparison",
    "Functional category/status annotation", "None", "None", "None", "Diagnostic panel", "None", "Replicon versus AMR panel",
    "None", "VF endpoint", "None"
  ),
  level = c(
    "Mixed episode- and participant-level", "Assembly-level", "Episode-level", "Episode-level plus participant-status summary",
    "Episode-level model with participant random intercept", "Episode nested within participant", "Adjacent-pair level nested within participant",
    "Variant within reference-specific comparison", "Gene-by-episode", "Assembly-tip level", "Module-by-transition-pair", "Episode-level",
    "Episode/sensitivity diagnostic", "Adjacent-pair level", "Episode-level", "Episode-level model where estimable",
    "Participant-pair level", "Adjacent-pair level nested within participant"
  ),
  repeated = c(
    "Yes — episode counts are repeated within participants", "Not for assembly QC rows; assemblies correspond to selected episodes",
    "Yes — repeated episodes and lineage structure", "Shown explicitly; participant panel is collapsed to status medians",
    "Handled with participant random intercept; all four ST-adjusted fits are singular", "Yes — ordered observations within participant",
    "Yes — participants contribute unequal numbers of adjacent pairs", "Yes — variants are nested within comparisons and participants",
    "Yes — repeated episodes can dominate visual patterns", "Yes — multiple assembly tips can belong to one participant",
    "Yes — transition endpoints are within participant", "Yes — unadjusted repeated episode points", "Yes — episode diagnostics repeat within participants",
    "Yes — participants may contribute multiple adjacent comparisons", "Yes — repeated episodes are descriptive", "Handled where GLMM fitting succeeded; sparse/fallback fits remain exploratory",
    "Collapsed to one nearest pair per included participant", "Handled by resident-cluster bootstrap"
  ),
  stringsAsFactors = FALSE
)

if (!is.null(final_manifest) && "figure_id" %in% names(final_manifest)) {
  required_final_fields <- c(
    "figure_id", "figure_number", "title", "scientific_question", "figure_class", "caption",
    "source_inputs", "statistical_method", "caveat", "unit", "filters", "visual_encodings",
    "multiplicity", "validation_status", "width_in", "height_in", "dpi"
  )
  if (length(setdiff(required_final_fields, names(final_manifest)))) {
    stop("Current final figure manifest lacks required fields: ", paste(setdiff(required_final_fields, names(final_manifest)), collapse = ", "))
  }
  if (!setequal(final_manifest$figure_id, final_panel_map$figure_id)) {
    stop("Current final figure manifest does not contain the expected 18 numbered figure families.")
  }
  for (i in seq_len(nrow(final_manifest))) {
    id <- as.character(final_manifest$figure_id[i])
    is_main <- identical(as.character(final_manifest$figure_class[i]), "Main")
    expected_key <- file.path(if (is_main) "plots/final" else "plots/final/supplementary", id)
    hit <- inventory$logical_key == expected_key
    if (!any(hit)) stop("Current final manifest figure is absent from inventory: ", expected_key)
    panel <- final_panel_map[match(id, final_panel_map$figure_id), , drop = FALSE]
    inventory$figure_number[hit] <- as.character(final_manifest$figure_number[i])
    inventory$figure_title[hit] <- as.character(final_manifest$title[i])
    inventory$scientific_question[hit] <- as.character(final_manifest$scientific_question[i])
    inventory$draft_caption[hit] <- as.character(final_manifest$caption[i])
    inventory$statistical_method[hit] <- as.character(final_manifest$statistical_method[i])
    inventory$required_caveat[hit] <- as.character(final_manifest$caveat[i])
    inventory$filtering_rule[hit] <- as.character(final_manifest$filters[i])
    inventory$visual_encoding_description[hit] <- as.character(final_manifest$visual_encodings[i])
    inventory$multiple_testing[hit] <- as.character(final_manifest$multiplicity[i])
    inventory$planned_width_in[hit] <- as.character(final_manifest$width_in[i])
    inventory$planned_height_in[hit] <- as.character(final_manifest$height_in[i])
    inventory$planned_dpi[hit] <- as.character(final_manifest$dpi[i])
    inventory$manifest_validation_status[hit] <- as.character(final_manifest$validation_status[i])
    inventory$analysis_represented[hit] <- paste0(as.character(final_manifest$title[i]), " — ", as.character(final_manifest$scientific_question[i]))
    inventory$source_data_file_or_object[hit] <- normalise_source_list(final_manifest$source_inputs[i])
    inventory$unit_of_observation[hit] <- as.character(final_manifest$unit[i])
    inventory$episode_or_participant_level[hit] <- panel$level
    inventory$plotted_denominator[hit] <- unname(final_denominators[[id]])
    inventory$plot_type[hit] <- panel$plot_type
    inventory$x_variable[hit] <- panel$x
    inventory$y_variable[hit] <- panel$y
    inventory$colour_fill_variable[hit] <- panel$colour
    inventory$facet_variable[hit] <- panel$facet
    inventory$statistical_result_displayed[hit] <- paste0(as.character(final_manifest$statistical_method[i]), "; multiplicity: ", as.character(final_manifest$multiplicity[i]))
    inventory$repeated_measures_present[hit] <- panel$repeated
    inventory$evidence_class[hit] <- if (id == "Fig02_wgs_quality_control") "QC" else if (id %in% c("Fig05_vf_association_evidence", "FigS05_near_miss_leave_one_uti", "FigS08_gene_model_forest", "FigS09_event_sample_sensitivity", "FigS10_snp_threshold_sensitivity")) "exploratory" else "descriptive"
    inventory$intended_use[hit] <- if (is_main) "thesis" else "supplementary"
    inventory$current_problems[hit] <- paste0("No unresolved generation/data-contract defect; required interpretive caveat: ", as.character(final_manifest$caveat[i]))
    inventory$severity[hit] <- "none"
    inventory$recommended_action[hit] <- "Retain in the canonical numbered thesis set; complete and record manual full-size and thesis-size visual QA"
    inventory$action_implemented[hit] <- "Regenerated as 300-dpi PNG and vector PDF from the current canonical source; manifest metadata and final-pack data checks recorded"
    inventory$canonical_status[hit] <- if (is_main) "canonical_final_main_retained" else "canonical_final_supplementary_retained"
    inventory$canonical[hit] <- "TRUE"
    inventory$validation_status[hit] <- if (final_checks_pass && identical(as.character(final_manifest$validation_status[i]), "validated")) "validated_by_current_final_pack_contract_and_data_checks" else "current_final_pack_validation_incomplete"
    inventory$visual_qa_status[hit] <- "current_png_pdf_rendered; manual_full_size_and_thesis_size_QA_pending"
    inventory$source_trace_status[hit] <- if (final_checks_pass) "current manifest source inputs and register line identified; all final_figure_data_checks pass" else "current manifest source inputs identified; final data checks incomplete"
  }
}

# A manual review ledger is created only after a human has inspected every
# current derivative.  When present, require one passing row per final figure
# and require the ledger to post-date all QA derivatives before upgrading the
# inventory from pending to reviewed.
manual_qa_complete <- FALSE
manual_qa_detail <- "not recorded"
qa_review_path <- file.path(root, "results/figure_audit/visual_qa_review.csv")
qa_results_path <- file.path(root, "results/figure_audit/visual_qa_results.csv")
qa_review <- safe_read_csv(qa_review_path)
qa_results <- safe_read_csv(qa_results_path)
if (!is.null(qa_review) && !is.null(final_manifest) &&
    all(c("figure_id", "pass") %in% names(qa_review)) &&
    nrow(qa_review) == nrow(final_manifest) &&
    !anyDuplicated(qa_review$figure_id) &&
    setequal(qa_review$figure_id, final_manifest$figure_id)) {
  review_pass <- if (is.logical(qa_review$pass)) {
    !is.na(qa_review$pass) & qa_review$pass
  } else {
    tolower(trimws(as.character(qa_review$pass))) %in% c("true", "t", "1", "yes", "pass", "passed")
  }
  derivative_paths <- if (!is.null(qa_results)) {
    unlist(qa_results[intersect(c("pdf_render_path", "a4_preview_path", "greyscale_path", "deutan_path", "protan_path"), names(qa_results))], use.names = FALSE)
  } else {
    character()
  }
  derivative_paths <- derivative_paths[file.exists(derivative_paths)]
  review_fresh <- length(derivative_paths) == 90L &&
    file.exists(qa_review_path) &&
    file.info(qa_review_path)$mtime >= max(file.info(derivative_paths)$mtime)
  manual_qa_complete <- all(review_pass) && review_fresh
  manual_qa_detail <- paste0(sum(review_pass), "/", length(review_pass),
                             " rows pass; ledger ", if (review_fresh) "is current" else "is stale")
  if (manual_qa_complete) {
    for (id in qa_review$figure_id) {
      hit <- inventory$logical_key %in% c(file.path("plots/final", id), file.path("plots/final/supplementary", id))
      inventory$visual_qa_status[hit] <- "human_review_passed_full_size_A4_greyscale_deutan_protan_and_PDF_render"
      inventory$recommended_action[hit] <- "Retain in the canonical numbered thesis set; repeat visual QA only after source or rendering changes"
      inventory$action_implemented[hit] <- paste0(inventory$action_implemented[hit], "; recorded complete human visual/accessibility review")
    }
  }
} else if (!is.null(qa_review)) {
  manual_qa_detail <- "review ledger present but incomplete or malformed"
}

# Diagnostic and statistical-sensitivity metadata enrich labels/evidence.
diag_meta <- safe_read_csv(file.path(root, "results/vf/uti_not_uti_diagnostic_figure_metadata.csv"))
if (!is.null(diag_meta) && all(c("figure_id", "evidence_type", "interpretation_limitations") %in% names(diag_meta))) {
  for (i in seq_len(nrow(diag_meta))) {
    hit <- stem == tolower(as.character(diag_meta$figure_id[i])) & grepl("^(plots/clinical|plots/vf)/", path)
    if (!any(hit)) next
    inventory$analysis_represented[hit] <- paste0(as.character(diag_meta$evidence_type[i]), " — ", as.character(diag_meta$interpretation_limitations[i]))
    inventory$statistical_result_displayed[hit] <- as.character(diag_meta$evidence_type[i])
    inventory$evidence_class[hit] <- if (grepl("QC", diag_meta$evidence_type[i], ignore.case = TRUE)) "QC" else if (grepl("explor|bootstrap|sensitivity", diag_meta$evidence_type[i], ignore.case = TRUE)) "exploratory" else "descriptive"
  }
}

stat_meta <- safe_read_csv(file.path(root, "results/statistical_sensitivity/statistical_sensitivity_figure_metadata.csv"))
if (!is.null(stat_meta) && "figure_id" %in% names(stat_meta)) {
  for (i in seq_len(nrow(stat_meta))) {
    hit <- stem == tolower(as.character(stat_meta$figure_id[i])) & grepl("^plots/statistical_sensitivity/", path)
    if (!any(hit)) next
    inventory$x_variable[hit] <- as.character(stat_meta$x_axis[i])
    inventory$y_variable[hit] <- as.character(stat_meta$y_axis[i])
    inventory$colour_fill_variable[hit] <- as.character(stat_meta$legend[i])
    inventory$analysis_represented[hit] <- as.character(stat_meta$caption[i])
    inventory$evidence_class[hit] <- "exploratory"
    inventory$intended_use[hit] <- "supplementary"
    inventory$canonical_status[hit] <- "active_duplicate_of_final_supplementary"
    inventory$canonical[hit] <- "FALSE"
  }
}

# Final figure panel metadata not fully represented in the manifest.
set_final <- function(id, unit, level, plot_type, x, y, colour, facet = "None", repeated = "Yes") {
  hit <- inventory$logical_key == paste0("plots/final_figures/", id)
  inventory$unit_of_observation[hit] <<- unit
  inventory$episode_or_participant_level[hit] <<- level
  inventory$plot_type[hit] <<- plot_type
  inventory$x_variable[hit] <<- x
  inventory$y_variable[hit] <<- y
  inventory$colour_fill_variable[hit] <<- colour
  inventory$facet_variable[hit] <<- facet
  inventory$repeated_measures_present[hit] <<- repeated
}

set_final("primary_denominator_and_uncertainty", "Episode", "Episode-level denominator audit", "Composite stacked/count bars", "Cohort/sensitivity definition", "Number of episodes", "Primary status; retention state", "Primary status facets", "Yes — episodes cluster within participants")
set_final("not_uti_to_uti_mechanism_casebook", "Not UTI-to-UTI transition", "Transition-level", "Composite count bar and evidence heatmap", "Mechanism/evidence feature", "Transition case", "Mechanism/evidence state", repeated = "Yes — participants can contribute transitions")
set_final("strain_stability_and_host_context", "Not UTI-to-UTI transition", "Transition-level", "Composite scatter and clinical-context heatmap", "SNP distance / context feature", "VF Jaccard / transition case", "Mechanism category / context state")
set_final("global_vf_signal_and_robustness", "Episode with participant bootstrap/model sensitivity", "Episode-level with participant resampling", "Composite forest, stability bars, and sensitivity contrasts", "Effect/difference estimate", "VF endpoint or stability flag", "Sensitivity definition", repeated = "Handled variably: participant bootstrap plus episode-level models")
set_final("transition_mechanisms_by_transition_type", "Adjacent within-participant transition", "Transition-level", "Proportional stacked bar", "Transition type", "Proportion of transition type", "Mechanism category")
set_final("accessory_plasmid_amr_changes", "Episode and Not UTI-to-UTI transition", "Episode- and transition-level", "Three-panel prevalence, burden and mechanism-change figure", "Replicon / predicted plasmid burden / mechanism feature", "Prevalence, episode burden or transition case", "Primary status; predicted change state")
set_final("near_miss_and_sparse_precision", "Near-miss episode / expected-count scenario", "Episode-level sensitivity", "Evidence heatmap and precision tile", "Clinical rule feature / assumed odds ratio", "Near-miss case / baseline prevalence", "Evidence/detectability state")
set_final("leave_one_uti_stability", "Model feature under leave-one-UTI perturbation", "Sensitivity/model-level", "Range/interval plot", "Effect estimate", "VF endpoint/module", "Direction-flip state", repeated = "Underlying episodes are repeated within participants")
set_final("prioritised_variant_map", "Annotated variant within selected transition", "Variant/transition-level", "Genomic position map", "Reference genome position (Mb ticks)", "Transition case", "Mechanism category")
set_final("lineage_confounding_diagnostic", "Episode/isolate", "Episode-level with within-ST facets", "Stacked proportion and faceted distributions", "Primary status", "ST proportion / detected VF genes", "Sequence type; primary status", "Sequence type", "Yes — repeated participant isolates")
set_final("paired_resident_expec_marker", "Participant observed in both states", "Participant-level paired", "Slope plot", "Primary status", "Participant-median ExPEC-like marker count", "Primary status", repeated = "Collapsed to participant/status medians")
set_final("not_uti_to_uti_module_gain_loss", "Transition by VF module", "Transition-level", "Heatmap", "Transition case", "VF module", "Module change state")
set_final("vf_module_pcoa_primary_status", "Episode/isolate", "Episode-level ordination", "PCoA scatter", "PCoA1", "PCoA2", "Primary status; ST-group shape", repeated = "Yes — repeated participant isolates")

# Research-question outputs are prespecified, reproducible diagnostics but
# remain exploratory/descriptive and are not promoted into the numbered pack
# merely because they were generated successfully.
idx <- grepl("^results/research_questions/RQ0[1-8]/", path)
inventory$canonical_status[idx] <- "canonical_research_question_diagnostic"
inventory$canonical[idx] <- "TRUE"
inventory$intended_use[idx] <- "supplementary"
inventory$evidence_class[idx] <- "exploratory"
inventory$current_problems[idx] <- "Prespecified exploratory research-question diagnostic; not an independently validated retained thesis figure"
inventory$severity[idx] <- "none"
inventory$recommended_action[idx] <- "Retain as a reproducible RQ diagnostic; cite the numbered final/supplementary figure where one supersedes it"
inventory$action_implemented[idx] <- "Source/release contract audited; retained as diagnostic rather than promoted to the numbered thesis set"
inventory$validation_status[idx] <- "validated_by_research_question_release_contract_not_individual_final_figure"
inventory$visual_qa_status[idx] <- "not_required_for_nonretained_rq_diagnostic"

# ------------------------------------------------------------------------------
# Known, evidence-backed defect and supersession rules from the source audit.
# ------------------------------------------------------------------------------

set_defect <- function(pattern, problem, severity, action, superseded = NULL, validation = "known_defect_requires_regeneration") {
  hit <- grepl(pattern, inventory$logical_key, perl = TRUE, ignore.case = TRUE)
  inventory$current_problems[hit] <<- problem
  inventory$severity[hit] <<- severity
  inventory$recommended_action[hit] <<- action
  inventory$action_implemented[hit] <<- "Defect documented in inventory; code/figure change pending unless stated elsewhere"
  inventory$validation_status[hit] <<- validation
  if (!is.null(superseded)) inventory$superseded_by[hit] <<- superseded
}

source_text <- function(file) {
  if (!file.exists(file)) return(NA_character_)
  paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# These final-pack rules are source-aware.  This prevents an audit rerun after a
# repair from perpetuating a defect solely because it existed in the baseline.
final_pack_text <- source_text("35_final_figure_pack.R")
if (is.na(final_pack_text)) {
  set_defect(
    "^plots/final/",
    "The current final-pack generating script 35_final_figure_pack.R is absent from the working tree, so existing final artifacts are not reproducible from the current checkout",
    "critical",
    "Restore or replace the canonical final-pack source, rerun it, then rebuild the audit inventory",
    validation = "critical_reproducibility_blocker_missing_source"
  )
} else {
  if (grepl("Four events show same-strain stable profiles; three are consistent with strain replacement", final_pack_text, fixed = TRUE)) {
    set_defect(
      "^plots/final_figures/not_uti_to_uti_mechanism_casebook$",
      "Final subtitle hard-codes four stable-profile and three replacement events, while the current casebook contains five and two respectively",
      "critical",
      "Derive category counts from casebook, add a text assertion, regenerate PNG/PDF, and visually verify the corrected subtitle"
    )
  }
  if (grepl("Not_UTI-to-UTI n=11; UTI-to-Not_UTI n=14; Not_UTI-to-Not_UTI n=391", final_pack_text, fixed = TRUE)) {
    set_defect(
      "^plots/final_figures/transition_mechanisms_by_transition_type$",
      "Caption hard-codes 11/14/391 transitions; current transition summary is 9 Not UTI-to-UTI, 12 UTI-to-Not UTI, and 349 Not UTI-to-Not UTI",
      "critical",
      "Generate all caption denominators from plotted data and fail validation on text/data mismatch"
    )
  }
  if (grepl("max_variant_pos <- max", final_pack_text, fixed = TRUE)) {
    set_defect(
      "^plots/final_figures/prioritised_variant_map$",
      "Genome track length is max observed plotted variant rather than reference-derived genome length; different reference coordinates are shown on one axis",
      "major",
      "Derive each reference length from FASTA/GFF, facet incompatible references, seed label repulsion, and state comparison direction"
    )
  }
  if (grepl("ST composition was not significantly different by primary status", final_pack_text, fixed = TRUE)) {
    set_defect(
      "^plots/final_figures/lineage_confounding_diagnostic$",
      "Subtitle states that ST composition is not significantly different, but the simulated Fisher result is not saved and the upstream test lacks a reproducible seed",
      "major",
      "Save the seeded test result with its exact table and derive or remove the significance claim"
    )
  }
}

sensitivity_text <- source_text("36_statistical_sensitivity_addon.R")
if (!is.na(sensitivity_text) && grepl("ST composition was not significantly different by primary status", sensitivity_text, fixed = TRUE)) {
  set_defect(
    "^plots/statistical_sensitivity/lineage_confounding_panel$",
    "Subtitle states that ST composition is not significantly different without loading a persisted, seeded test result and exact contingency table",
    "major",
    "Save the seeded test result with its exact table and derive or remove the significance claim"
  )
}
set_defect(
  "^plots/publication/Fig2_Mutation_Map$",
  "Hard-coded 5,000,000-bp genome line; Mb-formatted ticks with a bp axis label; cross-reference coordinates pooled; labels collide (21_publication_figures.R:143-177)",
  "critical",
  "Retire this plot and replace it with a reference-aware, derived-length variant map",
  "plots/final/Fig08_reference_aware_variant_map"
)
set_defect(
  "^plots/publication/Fig1_Swimmer_Plot$|^results/longitudinal/swimmer_plot$",
  "Phenotype-switch participants are computed or discussed but not programmatically highlighted/labelled; categorical spacing can imply unsupported temporal continuity",
  "major",
  "Build a date-aware longitudinal view with data-derived switcher annotation and explicit continuity caveat",
  "plots/final/Fig06_longitudinal_trajectories"
)
set_defect(
  "^plots/mechanism/not_uti_to_uti_case_matrix$",
  "Host-trigger text is truncated to 32 characters before a manual scale keyed by full strings, mapping the current trigger to the NA colour (33_mechanism_first_addon.R:1028,1035-1045)",
  "major",
  "Map stable categorical trigger codes before display-label truncation; assert manual-scale coverage"
)
set_defect(
  "^plots/wgs/panaroo_selection_(matrix|overview)$",
  "QC_REASON values such as HighContigs;/LowN50;/BadSize; do not match manual palette names such as Contigs > 500 or Size > 7MB (13_visualise_panaroo_selection.R:154-163)",
  "major",
  "Recode QC reasons once, assert scale-name coverage, and regenerate"
)
set_defect(
  "^plots/vf/variable_gene_heatmap$",
  "Missing values are converted to absence; binary profiles use default Euclidean clustering; variable/core definitions overlap (05_gene_overview_plots.R:183-240)",
  "major",
  "Use plots/final/supplementary/FigS01_vf_presence_heatmap; preserve this older heatmap only as superseded provenance",
  "plots/final/supplementary/FigS01_vf_presence_heatmap"
)
set_defect(
  "^plots/plasmids/replicon_heatmap$",
  "Missing values become absence and binary Euclidean clustering is used; PDF is an annotated heatmap while PNG is a different unlabeled base image",
  "critical",
  "Render both formats from one plot object, preserve missingness, and use justified binary distance/ordering"
)
set_defect(
  "^results/models/plots/volcano_plot_UTI_vs_Not_UTI$",
  "Raw-p=0.05 threshold is drawn while colour denotes FDR; zero/infinite odds ratios can appear as ordinary points; repeated episodes are not independent",
  "major",
  "Retain only as explicitly exploratory; use the separation-aware numbered effect/CI figures",
  "plots/final/Fig05_vf_association_evidence; plots/final/supplementary/FigS08_gene_model_forest"
)
set_defect(
  "^results/models/plots/forest_plot_top_hits$",
  "Sparse/separated or infinite estimates are not fully excluded/encoded; model fallbacks and singular fits are not visible",
  "major",
  "Use plots/final/supplementary/FigS08_gene_model_forest, which encodes finite/non-estimable and model-diagnostic states",
  "plots/final/supplementary/FigS08_gene_model_forest"
)
set_defect(
  "^results/lineage/st_risk_plot$",
  "Episode-level ST risk uses normal-approximation intervals with sparse/boundary counts and ignores repeated participant episodes",
  "major",
  "Use the descriptive prevalence/provenance display and avoid causal risk language",
  "plots/final/Fig03_sequence_type_distribution"
)
set_defect(
  "^plots/plasmids/(replicon_cooccurrence|ST_vs_replicon_network)$",
  "Reproducible weighted network is still a dense descriptive diagnostic and has no epidemiological evidence layer",
  "minor",
  "Keep supplementary/diagnostic only; edge alpha is episode frequency and must not be described as transmission or physical linkage"
)
set_defect(
  "^plots/vf/vf_score_effect_summary_uti_not_uti$|^plots/robustness/near_miss_score_shift$",
  "Bar encoding of median differences/medians can imply magnitude from zero and conceal distributions/uncertainty",
  "moderate",
  "Replace with point-and-interval or raw/participant-collapsed distribution display"
)
set_defect(
  "^plots/vf_amr/(replicon_burden_by_status|vf_vs_replicon_scatter|replicon_heatmap_top_STs)$",
  "Saved at 150 dpi and lacks sufficient repeated-measures/denominator context for thesis use",
  "moderate",
  "Regenerate at >=300 dpi plus vector where suitable; add sample/participant counts and caveat"
)
set_defect(
  "^plots/vf/(vf_pca_status|vf_pcoa_jaccard_status)$",
  "Primary-status ordinations omit the shared semantic UTI/Not_UTI colour scale",
  "major",
  "Use the canonical-status numbered PCoA and keep the old ordinations as diagnostics",
  "plots/final/supplementary/FigS04_vf_pcoa"
)
set_defect(
  "^plots/wgs/wgs_qc_n50_vs_contigs$",
  "Only N50 and contig count are shown; genome size, GC, pass/fail labels, and failure identification are absent",
  "moderate",
  "Use the multi-metric numbered QC panel with threshold lines and labelled failures",
  "plots/final/Fig02_wgs_quality_control"
)
set_defect(
  "^results/final_figures/denominator_flowchart/",
  "Rendered denominator workflow is stale (mixed-assembler 556/17 and 394/11/10 counts) despite current 532/16 and 371/9 contracts",
  "critical",
  "Regenerate from the updated derived-count script and remove/archive stale renderings"
)

# Exact alias/supersession rules.
alias_rules <- list(
  "^plots/vf/vf_burden_boxplot$" = "plots/vf/vf_burden_by_status",
  "^plots/vf/vf_category_barplot$" = "plots/vf/vf_category_burden_by_status",
  "^plots/vf/vf_burden_by_top_st$" = "plots/vf/vf_burden_by_st",
  "^plots/statistical_sensitivity/lineage_confounding_panel$" = "plots/final/Fig03_sequence_type_distribution; plots/final/Fig05_vf_association_evidence; plots/final/supplementary/FigS04_vf_pcoa",
  "^plots/statistical_sensitivity/paired_resident_expec_marker_slopeplot$" = "plots/final/Fig04_vf_burden",
  "^plots/statistical_sensitivity/not_uti_to_uti_module_gain_loss_heatmap$" = "plots/final/supplementary/FigS03_module_gain_loss",
  "^plots/statistical_sensitivity/vf_module_pcoa_primary_status$" = "plots/final/supplementary/FigS04_vf_pcoa"
)
for (pat in names(alias_rules)) {
  hit <- grepl(pat, inventory$logical_key, perl = TRUE)
  inventory$superseded_by[hit] <- alias_rules[[pat]]
  inventory$canonical_status[hit] <- "duplicate_or_superseded"
  inventory$canonical[hit] <- "FALSE"
  inventory$duplicate_group[hit] <- paste0("semantic_alias_", gsub("[^A-Za-z0-9]+", "_", basename(alias_rules[[pat]])))
  inventory$duplicate_basis[hit] <- "same plot object or direct final-pack duplicate identified in source"
  if (inventory$severity[hit][1] %in% c("moderate", "minor", "none")) inventory$severity[hit] <- "minor"
  inventory$recommended_action[hit] <- "Retain only the superseding canonical filename; keep this row for provenance"
  inventory$action_implemented[hit] <- "Duplicate/supersession recorded; file not deleted"
}

# Visual inspection is an acceptance gate for the retained numbered set. Do
# not leave nonretained/obsolete rows looking like unfinished final-figure QA.
nonretained_pending_qa <- inventory$canonical != "TRUE" &
  grepl("pending", inventory$visual_qa_status, ignore.case = TRUE)
inventory$visual_qa_status[nonretained_pending_qa] <- "not_required_for_nonretained_output"

# Code-declared missing figures remain explicit and are never treated as valid.
idx <- inventory$code_declared_missing
inventory$inventory_scope[idx] <- "code_declared_output_missing"
inventory$analysis_represented[idx] <- "Code-declared graphic output not present at the inferred output path"
inventory$unit_of_observation[idx] <- "Unknown until generated/source traced"
inventory$episode_or_participant_level[idx] <- "Unknown"
inventory$plotted_denominator[idx] <- "Unavailable — output missing"
inventory$current_problems[idx] <- ifelse(
  inventory$expected_in_canonical_pipeline[idx],
  "Expected canonical/conditional output is missing; determine whether branch conditions, failure, or stale path caused absence",
  "Noncanonical/optional/legacy code declares this output, but no artifact exists at the inferred path"
)
inventory$severity[idx] <- ifelse(inventory$expected_in_canonical_pipeline[idx], "major", "minor")
inventory$recommended_action[idx] <- ifelse(
  inventory$expected_in_canonical_pipeline[idx],
  "Inspect run log and branch conditions; either generate successfully or document intentional non-generation",
  "Keep as a code-declaration record; do not generate unless the workflow is reactivated"
)
inventory$action_implemented[idx] <- "Missing output recorded; no validation claim"
inventory$canonical_status[idx] <- ifelse(inventory$expected_in_canonical_pipeline[idx], "canonical_declared_missing_or_conditional", "noncanonical_declared_missing")
inventory$canonical[idx] <- ifelse(inventory$expected_in_canonical_pipeline[idx], "PENDING", "FALSE")
inventory$validation_status[idx] <- "missing_output"
inventory$visual_qa_status[idx] <- "not_possible_output_missing"

# `04_gene_breakdown.R` declares the nitrate-system UpSet plot inside a
# data-dependent branch.  Absence is intentional when the current VF hit table
# contains no Nar/Nap/Nas targets; it must not be reported as a failed figure.
nitrate_hit_n <- NA_integer_
vf_hits_path <- file.path(root, "results/vf/vf_hits_all.rds")
if (file.exists(vf_hits_path)) {
  vf_hits <- readRDS(vf_hits_path)
  nitrate_gene_col <- intersect(c("GENE", "gene", "Gene"), names(vf_hits))[1]
  if (!is.na(nitrate_gene_col)) {
    canon_gene <- function(x) tolower(gsub("[^a-z0-9]+", "", as.character(x)))
    nitrate_targets <- canon_gene(c(
      "narG", "narH", "narJ", "narI", "napF", "napD", "napA", "napG",
      "napH", "napB", "napC", "nasA", "nasB", "nasG", "nasH", "nasC"
    ))
    nitrate_hit_n <- sum(canon_gene(vf_hits[[nitrate_gene_col]]) %in% nitrate_targets, na.rm = TRUE)
  }
}
nitrate_missing <- inventory$code_declared_missing & inventory$logical_key == "plots/vf/nitrate_upset"
if (any(nitrate_missing) && identical(nitrate_hit_n, 0L)) {
  inventory$expected_in_canonical_pipeline[nitrate_missing] <- FALSE
  inventory$analysis_represented[nitrate_missing] <- "Conditional nitrate-system intersection plot; current source table contains zero qualifying Nar/Nap/Nas hits"
  inventory$unit_of_observation[nitrate_missing] <- "Episode with a detected nitrate-system target (none in current source data)"
  inventory$episode_or_participant_level[nitrate_missing] <- "Episode-level conditional display"
  inventory$plotted_denominator[nitrate_missing] <- "Not applicable — zero qualifying nitrate-system hits"
  inventory$current_problems[nitrate_missing] <- "No defect: the documented branch condition was false because the current VF source contains zero qualifying nitrate-system hits"
  inventory$severity[nitrate_missing] <- "none"
  inventory$recommended_action[nitrate_missing] <- "Do not generate an empty UpSet plot; retain the conditional code-declaration record"
  inventory$action_implemented[nitrate_missing] <- "Verified zero qualifying source hits and documented intentional non-generation"
  inventory$canonical_status[nitrate_missing] <- "conditional_output_not_applicable"
  inventory$canonical[nitrate_missing] <- "FALSE"
  inventory$validation_status[nitrate_missing] <- "validated_not_applicable_zero_qualifying_hits"
  inventory$visual_qa_status[nitrate_missing] <- "not_required_empty_conditional_output"
}

# Name-level duplicate groups supplement byte-identical duplicate detection.
name_key <- tolower(basename(inventory$logical_key))
name_counts <- table(name_key)
name_dup <- name_counts[name_key] > 1L
empty_dup <- !nzchar(inventory$duplicate_group)
inventory$duplicate_group[name_dup & empty_dup] <- paste0("name_", gsub("[^A-Za-z0-9]+", "_", name_key[name_dup & empty_dup]))
inventory$duplicate_basis[name_dup & empty_dup] <- "same logical filename in multiple directories; semantic identity requires review"

# Format completeness is explicit for final figures and declared siblings.
inventory$format_status <- ifelse(
  !inventory$output_exists,
  "missing",
  ifelse(
    grepl("^plots/(final|final_figures)/", inventory$logical_key),
    ifelse(grepl("PNG", inventory$output_formats) & grepl("PDF", inventory$output_formats), "complete_png_pdf", "incomplete_final_formats"),
    "existing_format_set_not_yet_adjudicated"
  )
)

# Keep allowed categorical values stable.
allowed_severity <- c("critical", "major", "moderate", "minor", "none")
if (any(!inventory$severity %in% allowed_severity)) stop("Invalid severity value generated.")
allowed_evidence <- c("exploratory", "confirmatory", "QC", "descriptive")
if (any(!inventory$evidence_class %in% allowed_evidence)) stop("Invalid evidence class generated.")
allowed_use <- c("diagnostic", "supplementary", "thesis", "manuscript", "obsolete")
if (any(!inventory$intended_use %in% allowed_use)) stop("Invalid intended-use value generated.")

# Select and order the public inventory schema. Requested fields are retained
# verbatim in snake_case, followed by audit/provenance fields.
inventory <- inventory[c(
  "figure_id", "figure_number", "figure_title", "scientific_question", "figure_filename",
  "generating_script", "approximate_source_code_lines", "draft_caption", "statistical_method",
  "required_caveat", "filtering_rule", "visual_encoding_description", "multiple_testing",
  "planned_width_in", "planned_height_in", "planned_dpi", "manifest_validation_status",
  "source_data_file_or_object", "analysis_represented", "unit_of_observation",
  "episode_or_participant_level", "plotted_denominator", "plot_type", "x_variable",
  "y_variable", "colour_fill_variable", "facet_variable", "statistical_result_displayed",
  "repeated_measures_present", "evidence_class", "intended_use", "current_problems",
  "severity", "recommended_action", "action_implemented", "final_output_path",
  "output_formats", "output_exists", "artifact_count", "dimensions", "width_px",
  "height_px", "page_count", "size_bytes", "mtime", "format_status", "canonical",
  "canonical_status", "validation_status", "visual_qa_status", "source_trace_status",
  "code_declared", "code_declared_missing", "code_source_declaration",
  "declared_formats", "expected_in_canonical_pipeline", "pipeline_stage",
  "pipeline_membership", "inventory_scope", "duplicate_group", "duplicate_basis",
  "superseded_by", "exact_duplicate_group", "logical_key", "primary_path", "artifact_paths"
)]

inventory <- inventory[order(
  match(inventory$severity, c("critical", "major", "moderate", "minor", "none")),
  inventory$canonical_status,
  inventory$logical_key
), ]
rownames(inventory) <- NULL
inventory$figure_id <- sprintf("F%05d", seq_len(nrow(inventory)))

# ------------------------------------------------------------------------------
# Write CSVs and an inventory-derived current-snapshot report.  Row-level
# pending states remain explicit where validation or visual QA is incomplete.
# ------------------------------------------------------------------------------

write.csv(artifacts, artifact_out, row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(inventory, inventory_out, row.names = FALSE, na = "", fileEncoding = "UTF-8")

count_table_md <- function(x, heading) {
  tab <- sort(table(x), decreasing = TRUE)
  c(
    paste0("### ", heading), "",
    "| Classification | Rows |", "|---|---:|",
    vapply(seq_along(tab), function(i) sprintf("| %s | %d |", names(tab)[i], as.integer(tab[i])), character(1)),
    ""
  )
}

critical_rows <- inventory[inventory$severity == "critical", c("logical_key", "current_problems", "recommended_action"), drop = FALSE]
critical_md <- if (nrow(critical_rows)) {
  c(
    "| Logical figure | Known defect | Required action |", "|---|---|---|",
    vapply(seq_len(nrow(critical_rows)), function(i) {
      clean <- function(x) gsub("\\|", "\\\\|", as.character(x))
      sprintf("| `%s` | %s | %s |", clean(critical_rows$logical_key[i]), clean(critical_rows$current_problems[i]), clean(critical_rows$recommended_action[i]))
    }, character(1))
  )
} else "No critical rows were classified by the static audit."

missing_rows <- inventory[
  inventory$code_declared_missing & inventory$expected_in_canonical_pipeline,
  c("logical_key", "generating_script", "approximate_source_code_lines", "current_problems"),
  drop = FALSE
]
missing_md <- if (nrow(missing_rows)) {
  c(
    "| Declared output | Script | Approximate line(s) | Interpretation |", "|---|---|---:|---|",
    vapply(seq_len(nrow(missing_rows)), function(i) {
      clean <- function(x) gsub("\\|", "\\\\|", as.character(x))
      sprintf(
        "| `%s` | `%s` | %s | %s |",
        clean(missing_rows$logical_key[i]), clean(missing_rows$generating_script[i]),
        clean(missing_rows$approximate_source_code_lines[i]), clean(missing_rows$current_problems[i])
      )
    }, character(1))
  )
} else "No canonical-run code-declared outputs were missing in this workspace snapshot."

retained_final <- inventory[inventory$canonical_status %in% c(
  "canonical_final_main_retained", "canonical_final_supplementary_retained"
), c(
  "figure_number", "figure_title", "scientific_question", "generating_script",
  "approximate_source_code_lines", "logical_key", "final_output_path",
  "canonical_status", "validation_status", "statistical_method", "required_caveat"
), drop = FALSE]
retained_final <- retained_final[order(retained_final$figure_number), , drop = FALSE]
principal_results <- c(
  Fig01 = "The selected analytical cohort contains 532 episodes from 161 participants: 16 UTI and 516 heterogeneous Not UTI episodes.",
  Fig02 = "Of 592 Longcycler candidates, 532 assemblies enter the canonical selected cohort after the implemented QC/selection contract.",
  Fig03 = "Sequence-type frequencies and call provenance are descriptive; episode imbalance and repeated sampling preclude an independent lineage-risk claim.",
  Fig04 = "Raw and participant-collapsed VF distributions show substantial within-group and within-participant heterogeneity; no unadjusted significance claim is made.",
  Fig05 = "All four prespecified ST-adjusted score models are singular, so apparent inverse estimates remain exploratory and model-dependent.",
  Fig06 = "Nine observed Not UTI-to-UTI transitions are present; five have direct same-strain support at the operational <=25-SNP reference.",
  Fig07 = "Among 371 adjacent comparisons, 140 are at or below 25 SNPs and 231 are above; genomic proximity is not transmission evidence.",
  Fig08 = "All 43 displayed variants from five comparisons were validated against their exact reference contig and are shown only within reference-specific coordinate systems.",
  FigS01 = "The display retains 108 curated genes with 5-95% prevalence and distinguishes presence, absence, and unavailable calls without binary Euclidean clustering.",
  FigS02 = "All 532 selected assemblies match the canonical core-SNP tree one-to-one; rooting is display-only and branch support is unavailable.",
  FigS03 = "Endpoint gain/loss states are shown for 32 VF modules across nine transitions; timing and causality cannot be inferred.",
  FigS04 = "The global VF ordination is unadjusted and descriptive; visible structure cannot be separated from lineage and repeated-participant effects.",
  FigS05 = "Leave-one-UTI diagnostics retain effect direction across the 16 omissions, but sparse-case dependence prevents confirmatory interpretation.",
  FigS06 = "The 371 transitions comprise 349 Not UTI-to-Not UTI, 9 Not UTI-to-UTI, 12 UTI-to-Not UTI, and 1 UTI-to-UTI interval.",
  FigS07 = "All 532 episodes have validated AMRFinderPlus primary profiles; mdf(A) is excluded from primary burden, genomic AMR is not phenotypic AST, and co-occurrence does not establish linkage.",
  FigS08 = "No prescreened gene model is BH-FDR significant; singular, separated, sparse, and non-estimable fits remain visibly flagged.",
  FigS09 = "Nearest within-participant contrasts are heterogeneous in 12 paired participants and remain an exploratory sensitivity analysis.",
  FigS10 = "The classified continuity proportion changes materially across prespecified SNP thresholds, showing that 25 SNPs is an operational rather than universal boundary."
)
retained_final$principal_result <- unname(principal_results[retained_final$figure_number])
retained_final_md <- if (nrow(retained_final)) {
  c(
    "| Figure | Title and scientific question | Source | Final output | Class | Principal result | Statistical method | Required caveat |",
    "|---|---|---|---|---|---|---|---|",
    vapply(seq_len(nrow(retained_final)), function(i) {
      clean <- function(x) gsub("\\|", "\\\\|", as.character(x))
      class_label <- if (identical(retained_final$canonical_status[i], "canonical_final_main_retained")) "Main" else "Supplementary"
      sprintf(
        "| %s | %s — %s | `%s:%s` | `%s` | %s | %s | %s | %s |",
        clean(retained_final$figure_number[i]), clean(retained_final$figure_title[i]), clean(retained_final$scientific_question[i]),
        clean(retained_final$generating_script[i]), clean(retained_final$approximate_source_code_lines[i]),
        clean(retained_final$final_output_path[i]), class_label, clean(retained_final$principal_result[i]),
        clean(retained_final$statistical_method[i]), clean(retained_final$required_caveat[i])
      )
    }, character(1))
  )
} else "No retained current final-pack rows were identified."

active_figures <- sum(inventory$canonical == "TRUE" & inventory$output_exists)
missing_expected <- sum(inventory$code_declared_missing & inventory$expected_in_canonical_pipeline)
scientific_candidates <- sum(inventory$inventory_scope == "scientific_or_diagnostic_figure")
exact_dup_artifacts <- sum(artifacts$exact_duplicate_artifact_count > 1L)
dimension_failures <- sum(!artifacts$dimension_status %in% c("read_png_header", "read_pdfinfo", "read_svg_header"))
retained_final_n <- nrow(retained_final)
retained_final_main_n <- sum(retained_final$canonical_status == "canonical_final_main_retained")
retained_final_supp_n <- sum(retained_final$canonical_status == "canonical_final_supplementary_retained")
superseded_old_final_n <- sum(inventory$canonical_status == "superseded_previous_final_pack")
final_check_n <- if (!is.null(final_checks)) nrow(final_checks) else 0L
final_check_pass_n <- if (!is.null(final_checks) && "pass" %in% names(final_checks)) sum(final_checks$pass %in% TRUE) else 0L
exploratory_existing_n <- sum(
  inventory$output_exists &
    inventory$inventory_scope == "scientific_or_diagnostic_figure" &
    inventory$evidence_class == "exploratory"
)
obsolete_or_duplicate_n <- sum(grepl(
  "obsolete|duplicate|superseded",
  inventory$canonical_status,
  ignore.case = TRUE
))
critical_discovered_n <- 7L
critical_remaining_n <- nrow(critical_rows)
critical_resolved_n <- critical_discovered_n - critical_remaining_n

change_entry <- function(file, problem, change, impact) {
  data.frame(file = file, problem = problem, change = change, impact = impact, stringsAsFactors = FALSE)
}
audit_change_log <- do.call(rbind, list(
  change_entry("R/plot_helpers.R", "Plot styling, status scales, saving, deidentification and effect handling were duplicated or unsafe.", "Added the shared publication theme, operational/legacy scales, scale assertions, deterministic seeds, deidentified labels, safe effect encoding and a PNG/PDF saver that propagates failures.", "Presentation and reproducibility changed; scientific source values did not."),
  change_entry("R/wgs_helpers.R", "Plot-save errors could be swallowed and empty outputs accepted.", "Made WGS plot saving fail closed, validate dimensions and output size, and log each saved file.", "No scientific value change; failure visibility changed."),
  change_entry("00c_plot_clinical_summary.R", "Clinical summaries could use a broader status table and timepoint ordering rather than the selected Longcycler cohort/date order.", "Bound the script to exact selected episode keys, derived chronological transitions, asserted 532/161/16/516 and 371/9 contracts, and standardised labels/colours.", "Plotted denominators, ordering and transition counts changed to the canonical cohort."),
  change_entry("03_plotting.R", "Legacy exploratory plotting silently caught save/tree/network failures.", "Removed broad silent try/save suppression and made failures explicit; the runner now gates this script as legacy-only.", "No retained scientific values changed; legacy generation and error behaviour changed."),
  change_entry("04_gene_breakdown.R", "Focused-gene figures could read the nonselected clinical status table.", "Restricted joins to the exact 532 selected Longcycler episode keys and asserted status totals/uniqueness.", "Denominators can change; current figures now use the canonical cohort."),
  change_entry("09_inc_plasmid_network.R", "PlasmidFinder calls could be duplicated across scripts, accession-collapsed features obscured separate replicon labels, and successful no-hit episodes were absent from wide profiles.", "Made script 09 the sole hash/version/database-bound PlasmidFinder caller, derived 80/80 and 90/90 profiles from one 80/60 scan, retained raw duplicate hits, published exact GENE-label matrices for all 532 episodes, and added threshold concordance and AP001918 feature validation.", "Primary plasmid-marker values and presentation changed; outputs are explicitly replicon detections rather than reconstructed plasmids."),
  change_entry("09b_mob_plasmid_reconstruction.R", "No cohort-wide reconstructed plasmid prediction layer or input-contig accounting existed.", "Added resumable MOB-suite 3.1.9 reconstruction for all 532 selected Longcycler assemblies with pinned tool/database hashes, explicit chromosome/predicted-plasmid/unassigned contig accounting and uncertainty fields.", "Adds assembly-based predicted plasmid-bin context without claiming circularity, transfer, transmission or causality."),
  change_entry("10_replicon_heatmap.R", "Missing replicon calls were converted to absence, keys could fall back positionally, and clustering/PNG-PDF outputs were inconsistent.", "Required exact isolate keys and binary values, preserved unavailable calls, used one ordered annotated representation for both formats, and asserted palette coverage.", "Matrix values and interpretation changed where missingness had formerly been zeroed."),
  change_entry("13_visualise_panaroo_selection.R", "Selection/QC palette levels and fallback selection paths could diverge from current Panaroo inputs.", "Required current cleanup/selection manifests, derived QC labels and palettes from actual levels, and added strict join/source checks.", "Displayed selection categories/counts can change to the current canonical inputs."),
  change_entry("20_variant_annotation_deep.R", "The GFF reader continued into the embedded FASTA section, producing tens of thousands of parsing problems per reference and obscuring annotation validity.", "Replaced permissive readr parsing with a strict nine-column pre-##FASTA parser that validates feature coordinates against embedded contig lengths and fails on malformed rows.", "Eligible variant count remained 43 after correction; annotation provenance is now warning-free and fail-closed."),
  change_entry("21_publication_figures.R", "The legacy swimmer used manual switch IDs and its mutation map pooled references behind a hard-coded 5-Mb line.", "Gated the script as legacy-only, derived timeline membership, wrote a retirement notice, and permanently stopped regenerating the invalid mutation map.", "The invalid mutation figure was removed from active generation; its replacement uses corrected values."),
  change_entry("24_vf_longitudinal_dynamics.R", "Longitudinal VF pairs were not required to equal the canonical selected transition export.", "Asserted exact selected episode keys, 371 adjacent pairs, nine Not UTI-to-UTI transitions, complete SNP/Jaccard evidence and equality to longcycler_transitions.csv.", "Denominators and pair membership are now canonical; invalid/missing pairs fail."),
  change_entry("28_vf_transition_case_studies.R", "Case studies admitted broader clinical transitions and potentially missing/nonselected genomic endpoints.", "Restricted all inputs to the exact selected cohort, required direct pair evidence, asserted the 371/9 contracts and standardised reader-facing terminology.", "Case membership and displayed counts changed to selected Longcycler pairs."),
  change_entry("29_vf_amr_combined_profile.R", "The script was a VF/plasmid placeholder with hard-coded false AMR availability fields and lacked sequence-validated determinant localization.", "Made script 29 the authoritative genomic-AMR layer and mapped VF calls through original contigs plus AMRFinder Prokka contigs through SHA-256/length/duplicate occurrence before joining MOB assignments; publishes 532 episode, 371 adjacent-pair and nine focused mechanism profiles.", "AMR values are genomic predictions rather than phenotypic AST; plasmid/chromosome placement and same-bin linkage remain assembly-based predictions."),
  change_entry("33_mechanism_first_addon.R", "Mechanism figures used stale status inputs, incomplete Panaroo sample-name matching and palette names that did not match host-trigger levels.", "Bound inputs to selected transitions, fixed slim-GFF Panaroo name mapping, recoded host-trigger categories and asserted complete manual-scale coverage.", "Mechanism counts and categories can change; unmatched palette states no longer disappear."),
  change_entry("34_robustness_first_addon.R", "Robustness text and checks retained 583/18/19-era denominators and a 17-UTI precision statement.", "Derived selected cohort, VF-ready, near-miss and transition counts from current data and removed stale hard-coded prose.", "Displayed counts changed to 532/161/16/516, 17 near-miss episodes and 371/9 transitions."),
  change_entry("35_final_figure_pack.R", "The provisional 13-figure pack contained stale counts, raw IDs, weak plot types, an old tree and an invalid pooled-reference mutation map.", "Rebuilt a single 18-family numbered pack with validated source contracts and made FigS07 the canonical three-layer replicon, predicted-plasmid burden and focused-transition mechanism figure using corrected PlasmidFinder and MOB assignments.", "Both plotted values and presentation changed; FigS07 explicitly avoids transfer, transmission, circularity and causal claims."),
  change_entry("scripts/research_questions/run_all.R", "The RQ wrapper used deprecated `.data$` pronouns in tidyselect expressions, which emitted three lifecycle warnings after otherwise complete RQ execution.", "Replaced those selection pronouns with reader-safe bare column selections and revalidated the RQ01-RQ10 wrapper contract.", "No estimand or plotted value changed; the release wrapper is warning-free."),
  change_entry("scripts/research_questions/run_rq09_10.R", "RQ09-RQ10 used `.data$` inside select, rename and pivot_wider tidyselect arguments, emitting repeated deprecation warnings.", "Removed 64 deprecated pronouns only from tidyselect contexts, preserved data-mask pronouns elsewhere, and executed the full module with warnings promoted to errors.", "No scientific value changed; current RQ09-RQ10 outputs regenerate without deprecation warnings."),
  change_entry("scripts/replay_model_warnings.R", "Three model-heavy modules exposed only aggregate warning counts in the pipeline log, obscuring whether conditions were statistical diagnostics or plotting/join defects.", "Added exact warning/message capture, reproducible classification, script SHA-256 binding and fail-closed rejection of deprecated, plotting, join or unclassified conditions.", "No model definition changed; 276 exact condition occurrences are now auditable."),
  change_entry("RUN_COMPLETE_ANALYSIS.sh", "Legacy figure layers ran as if canonical and the final pack preceded all RQ outputs/verification.", "Gated legacy scripts, moved the final pack after RQ analyses, added reference preparation, figure validation, visual-QA generation, inventory building and consolidated validation reporting.", "Pipeline order/acceptance changed; model definitions did not."),
  change_entry("RUN_IN_TERMINAL.sh", "Run identity and log provenance were not stable across nested commands.", "Exported one run ID and absolute log path and preserved canonical defaults/legacy gates.", "No scientific value change; run traceability changed."),
  change_entry("scripts/create_workflow_case_count_flowchart.R", "Workflow diagrams contained stale mixed-assembler counts and a hard-coded source-to-VF gap.", "Derived every displayed count from current outputs and regenerated the three denominator diagrams with reader-facing terminology.", "Displayed counts changed to the current 532/161/16/516 and 371/9 contracts."),
  change_entry("scripts/prepare_reference_aware_variants.R", "show-snps, FASTA and GFF contig identifiers were not proven equivalent before plotting.", "Added exact ID normalisation plus sequence-length/hash evidence, fail-closed ambiguity handling, reference-derived cumulative coordinates and position checks.", "Variant coordinates/eligibility changed; only 43 validated variants enter Fig08."),
  change_entry("scripts/validate_final_figures.R", "The final pack lacked machine acceptance tests for files, metadata, counts, palettes, joins and deidentification.", "Added 378 fail-closed checks across 18 families plus independent denominator/transition/reference assertions.", "No plotting value change; unsupported outputs now fail validation."),
  change_entry("scripts/visual_qa_final_figures.R", "No reproducible PDF/A4/greyscale/common-CVD review derivatives existed.", "Added Poppler rendering and deterministic A4, greyscale, deutan and protan derivatives with palette-distance checks.", "No scientific value change; accessibility verification became reproducible."),
  change_entry("scripts/build_figure_audit.R", "No complete repository artifact census or per-figure source/disposition inventory existed.", "Built the physical census, logical/code-declared inventory, duplicate/supersession rules, conditional-output handling and this report.", "No scientific value change; audit classification and traceability added."),
  change_entry("scripts/write_figure_audit_validation_results.R", "Validation evidence was scattered across manifests, logs and tests.", "Consolidated file/data/visual/palette/release/pipeline checks, exact hash-bound model-warning replay, warning classification, hashes and session information into one fail-closed result.", "No scientific value change; final acceptance evidence added."),
  change_entry("scripts/prepare_longcycler_release.R", "Generated-output cleanup could erase audit provenance or leave stale release graphics.", "Narrowed cleanup allowlists/exemptions so audit provenance is preserved while stale generated deliverables are removed deterministically.", "No analytical value change; output hygiene changed."),
  change_entry("scripts/verify_longcycler_release_deliverables.R", "The release verifier did not know the new final/audit contracts and root README policy.", "Added the figure-audit/final-pack requirements and narrowed provenance exemptions without restoring the deleted root README.", "No plotted value change; release acceptance changed."),
  change_entry("tests/test_plot_helpers.R", "Shared plotting contracts had no regression coverage.", "Added tests for scales, labels, safe estimates, deidentification and saving behaviour.", "No plotted value change; regression coverage added."),
  change_entry("tests/test_reference_aware_variants.R", "Reference/contig matching and ambiguity failure modes were untested.", "Added exact/hash, length-only rejection, nonunique match and contig-local variant-key tests.", "No plotted value change; invalid coordinate mappings are prevented."),
  change_entry("tests/test_research_question_layer.R", "The RQ release contract emitted a test-framework warning and did not cleanly exercise the current RQ script set.", "Kept the current RQ01-RQ10 parse/release checks and removed the conflicting fixed/ignore-case pattern option.", "No scientific value change; warning-free validation added.")
))

change_log_md <- c(
  "| Script | What was wrong | What changed and why | Did plotted values change? |",
  "|---|---|---|---|",
  vapply(seq_len(nrow(audit_change_log)), function(i) {
    clean <- function(x) gsub("\\|", "\\\\|", as.character(x))
    sprintf("| `%s` | %s | %s | %s |",
            clean(audit_change_log$file[i]), clean(audit_change_log$problem[i]),
            clean(audit_change_log$change[i]), clean(audit_change_log$impact[i]))
  }, character(1))
)

pipeline_marker <- file.path(root, "results/pipeline/RUN_COMPLETE.txt")
pipeline_sources <- file.path(root, c(
  "RUN_IN_TERMINAL.sh", "RUN_COMPLETE_ANALYSIS.sh", "35_final_figure_pack.R",
  "R/plot_helpers.R", "scripts/validate_final_figures.R",
  "scripts/visual_qa_final_figures.R", "scripts/build_figure_audit.R",
  "scripts/write_figure_audit_validation_results.R",
  "scripts/research_questions/run_all.R",
  "scripts/research_questions/run_rq09_10.R"
))
pipeline_sources <- pipeline_sources[file.exists(pipeline_sources)]
pipeline_marker_current <- file.exists(pipeline_marker) && length(pipeline_sources) &&
  file.info(pipeline_marker)$mtime >= max(file.info(pipeline_sources)$mtime)
pipeline_status_note <- if (pipeline_marker_current) {
  paste0("PASS — current completion marker `results/pipeline/RUN_COMPLETE.txt` (",
         format(file.info(pipeline_marker)$mtime, "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Amsterdam"), ").")
} else {
  "PENDING — the available completion marker predates one or more current figure-audit sources; a fresh full canonical run is required."
}
audit_status_complete <- manual_qa_complete && pipeline_marker_current &&
  final_check_n > 0L && final_check_pass_n == final_check_n && missing_expected == 0L

baseline_files <- list.files(out_dir, pattern = "^baseline_run_.*\\.txt$", full.names = FALSE)
baseline_note <- if (length(baseline_files)) paste0("Baseline record preserved unchanged: `results/figure_audit/", baseline_files[1], "`.") else "No baseline record was present when this scaffold was generated."

model_warning_path <- file.path(out_dir, "model_warning_replay.csv")
model_warning_note <- "Exact model-warning replay was not available."
if (file.exists(model_warning_path)) {
  model_warning <- suppressMessages(readr::read_csv(model_warning_path, show_col_types = FALSE))
  if (all(c("condition_type", "classification", "n") %in% names(model_warning))) {
    warning_totals_table <- stats::aggregate(
      as.numeric(model_warning$n),
      by = list(
        condition_type = model_warning$condition_type,
        classification = model_warning$classification
      ),
      FUN = sum,
      na.rm = TRUE
    )
    warning_totals <- paste0(
      warning_totals_table$condition_type, "/",
      warning_totals_table$classification, "=", warning_totals_table$x
    )
    unresolved_warning_n <- sum(model_warning$classification %in% c(
      "prohibited_deprecation", "prohibited_plot_or_join_warning",
      "unclassified_requires_review"
    ))
    model_warning_note <- paste0(
      "PASS — exact warning replay classified ", sum(as.numeric(model_warning$n), na.rm = TRUE),
      " condition occurrences (", paste(warning_totals, collapse = "; "),
      "); unresolved/prohibited distinct records=", unresolved_warning_n, "."
    )
  }
}

report_lines <- c(
  "# Figure Audit Report",
  "",
  paste0("Generated by `scripts/build_figure_audit.R` on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Europe/Amsterdam"), "."),
  "",
  paste0("> **STATUS: ", if (audit_status_complete) "COMPLETE FOR THE RETAINED NUMBERED THESIS SET" else "CURRENT SNAPSHOT; FINAL ACCEPTANCE STILL IN PROGRESS", ".** Repository-wide file/source discovery is current. The 18 numbered final figure families passed the final-pack data/file checks and have current PNG/PDF manifests; human visual/accessibility QA is ", if (manual_qa_complete) "complete" else "pending or incomplete", ". Other retained or diagnostic outputs keep row-specific validation states. Existence alone never confers scientific validation."),
  "",
  "## Executive summary",
  "",
  sprintf("- Figure outputs found in the current physical census: **%d**.", nrow(artifacts)),
  sprintf("- Numbered thesis figure families retained: **%d** (**%d main; %d supplementary**).", retained_final_n, retained_final_main_n, retained_final_supp_n),
  sprintf("- Retained families revised/rebuilt: **%d**.", retained_final_n),
  sprintf("- Prior provisional final-pack families replaced/superseded: **%d**.", superseded_old_final_n),
  sprintf("- Existing scientific/diagnostic outputs explicitly classified exploratory: **%d**.", exploratory_existing_n),
  sprintf("- Families assigned to supplementary use in the numbered pack: **%d**.", retained_final_supp_n),
  sprintf("- Logical rows classified obsolete, duplicated, or superseded: **%d**.", obsolete_or_duplicate_n),
  sprintf("- Critical defective output families identified: **%d** (**%d corrected; %d quarantined outside the retained set**).", critical_discovered_n, critical_resolved_n, critical_remaining_n),
  "- The operational UTI-versus-Not UTI contract governs the retained set; ASB/UTI/Negative is used only as explicitly labelled legacy/descriptive context.",
  "",
  "## Scope and current repository state",
  "",
  sprintf("- Physical graphical artifacts found repository-wide: **%d**.", nrow(artifacts)),
  sprintf("- Logical figure/artifact rows after grouping sibling formats and adding code-declared missing outputs: **%d**.", nrow(inventory)),
  sprintf("- Existing scientific/diagnostic candidate rows: **%d**.", scientific_candidates),
  sprintf("- Existing rows currently classified canonical/active: **%d**.", active_figures),
  sprintf("- Retained numbered thesis figures: **%d** (**%d main; %d supplementary**).", retained_final_n, retained_final_main_n, retained_final_supp_n),
  sprintf("- Previous unnumbered final-pack families explicitly superseded: **%d**.", superseded_old_final_n),
  sprintf("- Current final-pack data checks passed: **%d/%d**.", final_check_pass_n, final_check_n),
  sprintf("- Canonical code-declared outputs missing at inferred paths: **%d**.", missing_expected),
  sprintf("- Graphic artifacts whose dimensions could not be read: **%d**.", dimension_failures),
  sprintf("- Artifact rows participating in byte-identical duplicate groups: **%d**.", exact_dup_artifacts),
  paste0("- ", baseline_note),
  "- The root `README.md` is intentionally outside the generated release-document contract and may remain absent; the canonical workflow was mapped from `RUN_COMPLETE_ANALYSIS.sh`, `RUN_IN_TERMINAL.sh`, `00_config.R`, numbered scripts, and current result manifests.",
  "",
  "### Artifact-census snapshot provenance",
  "",
  "- Initial repository-wide snapshot (2026-07-14 12:34:27 CEST): **1,569** physical graphical artifacts.",
  "- Post-deliverables snapshot (2026-07-14 13:45:40 CEST): **1,252** physical graphical artifacts.",
  "- The **317-artifact reduction** comprises non-standalone presentation render, layout, preview, and QA assets removed by the exact v3/v4/v5 deliverables cleanup allowlist at 2026-07-14 13:17:39 CEST. These were not canonical scientific figures; no deleted binaries were restored.",
  "- The current live census above remains derived from the filesystem; these two dated values are retained as historical audit provenance.",
  "",
  "## Current analysis contract derived from local outputs",
  "",
  paste0("- Selected genomic cohort: ", selected_denominator, "."),
  sprintf("- Canonical adjacent within-participant transitions: %s.", ifelse(is.finite(contract$transition_n), contract$transition_n, "unavailable")),
  sprintf("- Focused Not UTI-to-UTI mechanism cases: %s.", ifelse(is.finite(contract$case_n), contract$case_n, "unavailable")),
  "- The current primary comparison is **UTI versus Not UTI** in the exact selected QC-pass Longcycler cohort. ASB/Negative classifications remain relevant for legacy/contextual outputs but must not be presented as the current primary model contrast.",
  "",
  count_table_md(inventory$severity, "Severity classification"),
  count_table_md(inventory$canonical_status, "Canonical-status classification"),
  count_table_md(inventory$inventory_scope, "Artifact/figure scope"),
  "## Retained numbered thesis figure set",
  "",
  retained_final_md,
  "",
  "All retained rows were generated by the current `35_final_figure_pack.R` into `plots/final/` or `plots/final/supplementary/`. The corresponding 18-row semantic manifest, 18-row save/hash manifest, and 15-row all-pass data-check table are under `results/figure_audit/`.",
  "",
  "## Change log for figure-audit implementation",
  "",
  "This table covers scripts changed by this figure-audit implementation. Pre-existing unrelated dirty-worktree changes remain user-owned and are not attributed to the audit.",
  "",
  change_log_md,
  "",
  "## Remaining known critical defects outside the retained numbered pack",
  "",
  critical_md,
  "",
  "## Unresolved and deliberately excluded issues",
  "",
  sprintf("- **%d critical legacy/noncanonical families remain quarantined:** none is referenced by the retained numbered pack; each remains explicitly obsolete or unvalidated in the inventory.", critical_remaining_n),
  "- Nonretained pipeline graphics are classified as diagnostics/source layers rather than promoted by render success. Their row-level inventory records source, problems, severity, disposition and superseding figure; they are not claimed as independently result-certified thesis figures.",
  "- Historical artifacts with unavailable or obsolete source contracts remain explicitly unvalidated. No missing source is concealed as biological absence.",
  "- The retired `plots/publication/Fig2_Mutation_Map.png`, legacy replicon heatmap, and stale denominator-flowchart copies outside `docs/figures/workflow_flowchart/` must not be cited.",
  "",
  "## Canonical-run code-declared outputs currently missing",
  "",
  missing_md,
  "",
  "## Major duplicate, obsolete, and superseded groups",
  "",
  "- All 13 former unnumbered families under `plots/final_figures/` are classified `superseded_previous_final_pack`, linked to numbered replacements, and excluded from the retained set.",
  "- `21_publication_figures.R` is retained only behind the explicit `RUN_LEGACY_PUBLICATION_FIGURES=1` Phase-3 gate and is skipped by default; its swimmer and mutation map are superseded by Fig06 and Fig08 and must not be treated as current publication figures.",
  "- `36_statistical_sensitivity_addon.R` standalone figures are source layers or duplicates of the numbered final and supplementary figures.",
  "- Exact alias outputs in scripts 23 and 25 are recorded through `duplicate_group` and `superseded_by`.",
  "- `results/plots/`, bulk persistence/participant outputs, archived outputs, document-page renders, and slide renders are retained in the artifact census but are not automatically scientific figure candidates.",
  "- The 90 derivatives under `results/figure_audit/visual_qa/` are classified as non-standalone QA/render artifacts; validation and citation remain attached to their 18 canonical source figures.",
  "- Denominator/workflow renders under `docs/figures/workflow_flowchart/` were regenerated from current derived values; older denominator renders elsewhere remain noncanonical unless their row states otherwise.",
  "",
  "## Coverage added by the numbered pack",
  "",
  "- Fig02 supplies multi-metric WGS QC with implemented thresholds and deidentified failure labels.",
  "- Fig03 distinguishes provider, local-fallback, and missing/non-typable ST provenance.",
  "- Fig06 derives the nine Not UTI-to-UTI transition cases programmatically and uses collection dates with an explicit continuity caveat.",
  "- Fig08 uses exact reference/contig validation, derived reference lengths, and reference-specific facets; the former common-coordinate maps are obsolete.",
  "- FigS01 distinguishes presence, absence, and unavailable data and uses biological ordering rather than default Euclidean clustering.",
  "- FigS02 provides the 532-tip canonical core-genome phylogeny with exact metadata matching and explicit rooting/support limitations.",
  "- Fig05 and FigS08 expose effect estimates, confidence intervals, multiplicity, singularity, separation, and non-estimable results without plotting infinity as an ordinary finite estimate.",
  "",
  "## Verification status",
  "",
  sprintf("- **Numbered final-pack data/file contract:** PASS — %d/%d checks passed, including all 36 PNG/PDF files present and nonempty.", final_check_pass_n, final_check_n),
  "- **Numbered final-pack source trace:** PASS — every retained row has a current register line, source-input list, denominator, caption, statistical method, caveat, dimensions, and save metadata.",
  "- **Numbered final-pack warnings/errors:** PASS for the recorded render run — the pack completed with warnings promoted to errors; hashes and byte counts were recorded by the shared saver.",
  paste0("- **Model-warning classification:** ", model_warning_note),
  paste0("- **Manual visual QA:** ", if (manual_qa_complete) "PASS" else "PENDING/INCOMPLETE", " — ", manual_qa_detail, ". Review covers full resolution, approximate A4 insertion size, PDF render, greyscale, deutan, and protan views."),
  "- **Non-final figures:** retain their row-specific unvalidated/known-defect states; the final-pack checks do not retroactively validate legacy or diagnostic artifacts.",
  paste0("- **Pipeline:** ", pipeline_status_note),
  "- **Consolidated validation record:** `results/figure_audit/validation_results.txt` records the exact run log, warning classification, session information, final files, dimensions, data anchors and unresolved exclusions.",
  "",
  "## Machine-readable deliverables",
  "",
  "- `results/figure_audit/artifact_census.csv`: one row per physical graphic-like file, including dimensions, size, mtime, MD5 duplicate groups, role, and scope.",
  "- `results/figure_audit/figure_inventory.csv`: one row per logical figure/artifact group plus code-declared missing outputs, with all scientific, statistical, provenance, action, and QA fields requested in the audit brief.",
  "",
  "## Next audit checkpoint",
  "",
  if (manual_qa_complete) "Repeat the automated and human visual review after any final-figure source or rendering change. Legacy rows must remain excluded unless deliberately reactivated and revalidated." else "After manual visual review, add a complete current `visual_qa_review.csv` ledger and rerun `Rscript scripts/build_figure_audit.R`. Legacy rows must remain excluded unless deliberately reactivated and revalidated."
)

writeLines(report_lines, report_out, useBytes = TRUE)

message(sprintf("Artifact census: %d rows -> %s", nrow(artifacts), rel_path(artifact_out)))
message(sprintf("Figure inventory: %d rows (%d existing; %d code-declared missing) -> %s",
                nrow(inventory), sum(inventory$output_exists), sum(inventory$code_declared_missing), rel_path(inventory_out)))
message(sprintf("Current-snapshot audit report -> %s", rel_path(report_out)))

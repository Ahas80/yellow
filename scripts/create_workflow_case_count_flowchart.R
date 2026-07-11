#!/usr/bin/env Rscript

# Build the numbered workflow-denominator figures and Markdown explainer.
# The output is documentation, not a pipeline input.

root <- normalizePath(".", winslash = "/", mustWork = TRUE)

paths <- list(
  cohort_flow = file.path(root, "results/summary/table_01_cohort_episode_flow.csv"),
  denominator = file.path(root, "results/qc/pipeline_denominator_summary.csv"),
  transition_index = file.path(root, "results/vf/vf_transition_case_index.csv"),
  vf_longitudinal = file.path(root, "results/vf/vf_longitudinal_transitions.csv"),
  out_dir = file.path(root, "docs/figures/workflow_flowchart"),
  md = file.path(root, "docs/workflow_case_count_flowchart.md")
)

for (path in paths[c("cohort_flow", "denominator", "transition_index", "vf_longitudinal")]) {
  if (!file.exists(path)) {
    stop("Required source file is missing: ", path, call. = FALSE)
  }
}

dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

cohort_flow <- read_csv(paths$cohort_flow)
denominator <- read_csv(paths$denominator)
transition_index <- read_csv(paths$transition_index)
vf_longitudinal <- read_csv(paths$vf_longitudinal)

pick_denominator <- function(stage) {
  row <- denominator[denominator$stage == stage, , drop = FALSE]
  if (nrow(row) == 0) {
    stop("Missing denominator stage: ", stage, call. = FALSE)
  }
  row[nrow(row), , drop = FALSE]
}

pick_flow <- function(stage_order) {
  row <- cohort_flow[cohort_flow$stage_order == stage_order, , drop = FALSE]
  if (nrow(row) == 0) {
    stop("Missing cohort-flow stage_order: ", stage_order, call. = FALSE)
  }
  row[1, , drop = FALSE]
}

truthy <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

fmt_n <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

status_all <- pick_denominator("status_map")
status_included <- pick_denominator("status_map_primary_included")
assembly_qc <- pick_denominator("assembly_qc")
canonical <- pick_denominator("canonical_selected_assemblies")
vf_ready <- pick_denominator("vf_analysis_ready")
model_dataset <- pick_denominator("model_dataset")
longitudinal_subset <- pick_flow(6)

transition_breakdown <- table(transition_index$transition_type, useNA = "no")
clinical_not_uti_to_uti <- sum(truthy(transition_index$is_not_uti_to_uti))
linked_not_uti_to_uti <- sum(
  truthy(transition_index$is_not_uti_to_uti) & truthy(transition_index$has_vf_pair)
)
missing_linked_not_uti_to_uti <- clinical_not_uti_to_uti - linked_not_uti_to_uti

counts <- list(
  classified_episodes = as.integer(status_all$n_rows),
  classified_participants = as.integer(status_all$n_unique_participants),
  classified_uti = as.integer(status_all$n_UTI),
  classified_not_uti = as.integer(status_all$n_Not_UTI),
  included_episodes = as.integer(status_included$n_rows),
  included_participants = as.integer(status_included$n_unique_participants),
  included_uti = as.integer(status_included$n_UTI),
  included_not_uti = as.integer(status_included$n_Not_UTI),
  assembly_records = as.integer(assembly_qc$n_rows),
  canonical_episodes = as.integer(canonical$n_rows),
  canonical_participants = as.integer(canonical$n_unique_participants),
  vf_ready_episodes = as.integer(vf_ready$n_rows),
  vf_ready_participants = as.integer(vf_ready$n_unique_participants),
  vf_ready_uti = as.integer(vf_ready$n_UTI),
  vf_ready_not_uti = as.integer(vf_ready$n_Not_UTI),
  model_episodes = as.integer(model_dataset$n_rows),
  longitudinal_episodes = as.integer(longitudinal_subset$n_episodes),
  longitudinal_participants = as.integer(longitudinal_subset$n_participants),
  longitudinal_comparisons = nrow(vf_longitudinal),
  longitudinal_comparison_participants = length(unique(vf_longitudinal$Participant_id)),
  total_transitions = nrow(transition_index),
  transition_participants = length(unique(transition_index$Participant_id)),
  clinical_not_uti_to_uti = clinical_not_uti_to_uti,
  linked_not_uti_to_uti = linked_not_uti_to_uti,
  missing_linked_not_uti_to_uti = missing_linked_not_uti_to_uti
)

checks <- data.frame(
  item = c(
    "included clinical visits",
    "included participants",
    "included UTI",
    "included Not_UTI",
    "assembly candidates",
    "selected genomic profiles",
    "selected genomic participants",
    "VF-ready profiles",
    "VF-ready UTI",
    "VF-ready Not_UTI",
    "longitudinal visit subset",
    "longitudinal subset participants",
    "longitudinal comparisons",
    "status comparisons",
    "clinical Not_UTI -> UTI transitions",
    "WGS/VF-linked Not_UTI -> UTI transitions"
  ),
  observed = c(
    counts$included_episodes,
    counts$included_participants,
    counts$included_uti,
    counts$included_not_uti,
    counts$assembly_records,
    counts$canonical_episodes,
    counts$canonical_participants,
    counts$vf_ready_episodes,
    counts$vf_ready_uti,
    counts$vf_ready_not_uti,
    counts$longitudinal_episodes,
    counts$longitudinal_participants,
    counts$longitudinal_comparisons,
    counts$total_transitions,
    counts$clinical_not_uti_to_uti,
    counts$linked_not_uti_to_uti
  ),
  expected = c(583, 166, 18, 565, 1291, 556, 162, 556, 17, 539, 538, 144, 394, 417, 11, 10)
)

bad_checks <- checks[checks$observed != checks$expected, , drop = FALSE]
if (nrow(bad_checks) > 0) {
  msg <- apply(
    bad_checks,
    1,
    function(row) paste0(row[["item"]], ": observed ", row[["observed"]], ", expected ", row[["expected"]])
  )
  stop(paste(c("Count checks failed:", msg), collapse = "\n"), call. = FALSE)
}

library(grid)

palette <- list(
  ink = "#0F172A",
  muted = "#475569",
  faint = "#64748B",
  line = "#CBD5E1",
  pale = "#F8FAFC",
  blue = "#0072B2",
  blue_fill = "#E0F2FE",
  orange = "#D55E00",
  orange_fill = "#FFF7ED",
  green = "#009E73",
  green_fill = "#ECFDF5",
  slate_fill = "#F1F5F9",
  yellow = "#E69F00",
  yellow_fill = "#FEF9C3",
  white = "#FFFFFF"
)

wrap_text <- function(text, width) {
  paste(strwrap(text, width = width), collapse = "\n")
}

text_grob <- function(label, x, y, width, size = 10, color = palette$ink,
                      bold = FALSE, just = c("left", "top")) {
  grid.text(
    label,
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    just = just,
    gp = gpar(
      fontsize = size,
      col = color,
      fontface = if (bold) "bold" else "plain",
      lineheight = 0.98
    )
  )
}

draw_title <- function(title, subtitle) {
  grid.rect(gp = gpar(fill = palette$white, col = NA))
  text_grob(title, 0.045, 0.955, 0.9, size = 21, bold = TRUE)
  text_grob(subtitle, 0.045, 0.905, 0.9, size = 10.5, color = palette$faint)
  grid.lines(
    x = unit(c(0.045, 0.955), "npc"),
    y = unit(c(0.855, 0.855), "npc"),
    gp = gpar(col = palette$line, lwd = 1)
  )
}

draw_card <- function(x, y, w, h, idx, title, metric, detail,
                      fill = palette$slate_fill, line = palette$line,
                      accent = palette$ink) {
  grid.roundrect(
    x = unit(x + w / 2, "npc"),
    y = unit(y - h / 2, "npc"),
    width = unit(w, "npc"),
    height = unit(h, "npc"),
    r = unit(0.018, "npc"),
    gp = gpar(fill = fill, col = line, lwd = 1.1)
  )
  grid.roundrect(
    x = unit(x + 0.024, "npc"),
    y = unit(y - 0.036, "npc"),
    width = unit(0.036, "npc"),
    height = unit(0.046, "npc"),
    r = unit(0.006, "npc"),
    gp = gpar(fill = accent, col = NA)
  )
  grid.text(
    idx,
    x = unit(x + 0.024, "npc"),
    y = unit(y - 0.036, "npc"),
    gp = gpar(fontsize = 8.3, col = palette$white, fontface = "bold")
  )
  text_grob(wrap_text(title, 28), x + 0.050, y - 0.016, w - 0.066, size = 8.8, bold = TRUE)
  text_grob(metric, x + 0.020, y - 0.075, w - 0.040, size = 17, color = accent, bold = TRUE)
  text_grob(wrap_text(detail, 30), x + 0.020, y - 0.122, w - 0.040, size = 7.7, color = palette$faint)
}

draw_arrow <- function(x1, y1, x2, y2, color = "#94A3B8") {
  grid.segments(
    x0 = unit(x1, "npc"),
    y0 = unit(y1, "npc"),
    x1 = unit(x2, "npc"),
    y1 = unit(y2, "npc"),
    gp = gpar(col = color, lwd = 1.6),
    arrow = grid::arrow(type = "closed", length = unit(0.014, "npc"))
  )
}

draw_note <- function(x, y, label, color = palette$faint) {
  text_grob(wrap_text(label, 42), x, y, 0.28, size = 7.3, color = color)
}

draw_label_box <- function(label, x, y, w, h = 0.040) {
  grid.roundrect(
    x = unit(x + w / 2, "npc"),
    y = unit(y - h / 2, "npc"),
    width = unit(w, "npc"),
    height = unit(h, "npc"),
    r = unit(0.010, "npc"),
    gp = gpar(fill = palette$white, col = palette$line, lwd = 0.7)
  )
  text_grob(wrap_text(label, 42), x + 0.010, y - 0.010, w - 0.020, size = 6.9, color = palette$muted)
}

draw_unit_story <- function() {
  draw_title(
    "Why the counts change: each stage counts a different thing",
    "Clinical visits expand into assembly candidates for QC, then collapse back to one selected genomic profile per usable visit."
  )

  text_grob("1. Clinical curation", 0.055, 0.742, 0.18, size = 10, color = palette$orange, bold = TRUE)
  draw_card(0.230, 0.770, 0.220, 0.135, "1", "Classified clinical visits",
            paste0(fmt_n(counts$classified_episodes), " visits"),
            "participant sample visits",
            fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange)
  draw_card(0.555, 0.770, 0.240, 0.135, "2", "Included clinical visits",
            paste0(fmt_n(counts$included_episodes), " visits"),
            paste0(fmt_n(counts$included_uti), " UTI; ",
                   fmt_n(counts$included_not_uti), " Not_UTI"),
            fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange)

  draw_arrow(0.455, 0.705, 0.548, 0.705)
  draw_label_box("-2 manual exclusions", 0.455, 0.755, 0.110)

  text_grob("2. Temporary assembly inventory", 0.055, 0.552, 0.20, size = 10, color = palette$blue, bold = TRUE)
  draw_card(0.370, 0.575, 0.300, 0.135, "3", "Assembly candidates / QC",
            paste0(fmt_n(counts$assembly_records), " candidates"),
            "genome files, not extra visits",
            fill = palette$blue_fill, line = "#BAE6FD", accent = palette$blue)
  draw_arrow(0.675, 0.635, 0.675, 0.585)
  draw_label_box("counting changes: visits -> assembly candidates", 0.695, 0.642, 0.220)

  text_grob("3. Return to one profile per visit", 0.055, 0.382, 0.22, size = 10, color = palette$green, bold = TRUE)
  draw_card(0.370, 0.405, 0.300, 0.135, "4", "Selected genomic profiles",
            paste0(fmt_n(counts$canonical_episodes), " profiles"),
            "one profile per usable visit",
            fill = palette$green_fill, line = "#BBF7D0", accent = palette$green)
  draw_arrow(0.520, 0.440, 0.520, 0.405)
  draw_label_box("QC selects one profile per usable visit", 0.695, 0.455, 0.230)

  text_grob("4. Analysis-specific layers", 0.055, 0.220, 0.20, size = 10, color = palette$muted, bold = TRUE)
  draw_card(0.215, 0.250, 0.220, 0.115, "5A", "VF/model-ready table",
            paste0(fmt_n(counts$vf_ready_episodes), " profiles"),
            paste0(fmt_n(counts$vf_ready_uti), " UTI; ",
                   fmt_n(counts$vf_ready_not_uti), " Not_UTI"),
            fill = palette$slate_fill, line = palette$line, accent = palette$green)
  draw_card(0.475, 0.250, 0.220, 0.115, "5B", "Longitudinal visit layer",
            paste0(fmt_n(counts$longitudinal_episodes), " visits"),
            paste0(fmt_n(counts$longitudinal_comparisons), " consecutive visit pairs"),
            fill = palette$slate_fill, line = palette$line, accent = palette$blue)
  draw_card(0.735, 0.250, 0.220, 0.115, "5C", "Focused transition layer",
            paste0(fmt_n(counts$clinical_not_uti_to_uti), " transitions"),
            paste0(fmt_n(counts$linked_not_uti_to_uti), " WGS/VF-linked"),
            fill = palette$slate_fill, line = palette$line, accent = palette$orange)

  draw_arrow(0.520, 0.270, 0.325, 0.252)
  draw_arrow(0.520, 0.270, 0.585, 0.252)
  draw_arrow(0.520, 0.270, 0.845, 0.252)

  grid.roundrect(
    x = unit(0.500, "npc"),
    y = unit(0.065, "npc"),
    width = unit(0.865, "npc"),
    height = unit(0.060, "npc"),
    r = unit(0.016, "npc"),
    gp = gpar(fill = palette$pale, col = palette$line, lwd = 0.9)
  )
  text_grob(
    wrap_text("Mental model: 585/583 are clinical visits; 1,291 is a temporary list of assembly candidates; 556 is selected genomic profiles; 538, 394, 11, and 10 are downstream analysis-specific subsets.", 148),
    0.085, 0.083, 0.80, size = 7.8, color = palette$muted
  )
  text_grob("Sources: pipeline_denominator_summary.csv; table_01_cohort_episode_flow.csv; vf_longitudinal_transitions.csv; vf_transition_case_index.csv",
            0.045, 0.018, 0.9, size = 6.9, color = palette$faint)
}

draw_ladder <- function() {
  draw_title(
    "Numbered workflow denominators with plain-language units",
    "Each card title names the unit; the large number shows the count for that stage."
  )

  w <- 0.168
  h <- 0.194
  y <- 0.760
  xs <- c(0.040, 0.232, 0.424, 0.616, 0.808)
  centers_y <- y - h / 2

  draw_card(xs[1], y, w, h, "1", "Classified clinical visits",
            paste0(fmt_n(counts$classified_episodes), " visits"),
            paste0(fmt_n(counts$classified_participants), " participants; ",
                   fmt_n(counts$classified_uti), " UTI, ",
                   fmt_n(counts$classified_not_uti), " Not_UTI before exclusions"),
            fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange)
  draw_card(xs[2], y, w, h, "2", "Included clinical visits",
            paste0(fmt_n(counts$included_episodes), " visits"),
            paste0(fmt_n(counts$included_participants), " participants; ",
                   fmt_n(counts$included_uti), " UTI, ",
                   fmt_n(counts$included_not_uti), " Not_UTI; 2 exclusions"),
            fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange)
  draw_card(xs[3], y, w, h, "3", "Assembly candidates for QC",
            fmt_n(counts$assembly_records),
            "Assembler alternatives stay explicit here; these are not extra visits",
            fill = palette$blue_fill, line = "#BAE6FD", accent = palette$blue)
  draw_card(xs[4], y, w, h, "4", "Selected genomic profiles",
            paste0(fmt_n(counts$canonical_episodes), " profiles"),
            paste0(fmt_n(counts$canonical_participants),
                   " participants; VF calls; one profile per usable visit"),
            fill = palette$blue_fill, line = "#BAE6FD", accent = palette$blue)
  draw_card(xs[5], y, w, h, "5", "VF/model-ready profiles",
            paste0(fmt_n(counts$vf_ready_episodes), " profiles"),
            paste0(fmt_n(counts$vf_ready_uti), " UTI, ",
                   fmt_n(counts$vf_ready_not_uti),
                   " Not_UTI; 27 included visits lack VF-ready evidence"),
            fill = palette$green_fill, line = "#BBF7D0", accent = palette$green)

  for (i in 1:4) {
    draw_arrow(xs[i] + w + 0.003, centers_y, xs[i + 1] - 0.003, centers_y)
  }

  branch_y <- 0.335
  draw_card(0.222, branch_y, 0.256, 0.174, "6", "Longitudinal visit subset",
            paste0(fmt_n(counts$longitudinal_episodes), " visits"),
            paste0(fmt_n(counts$longitudinal_participants), " participants; ",
                   fmt_n(counts$longitudinal_comparisons), " consecutive visit pairs"),
            fill = palette$slate_fill, line = palette$line, accent = palette$blue)
  draw_card(0.544, branch_y, 0.256, 0.174, "7", "Not_UTI -> UTI transitions",
            paste0(fmt_n(counts$clinical_not_uti_to_uti), " transitions"),
            paste0(fmt_n(counts$linked_not_uti_to_uti),
                   " WGS/VF-linked; ",
                   fmt_n(counts$missing_linked_not_uti_to_uti), " missing endpoint"),
            fill = palette$slate_fill, line = palette$line, accent = palette$orange)

  draw_arrow(xs[5] + w / 2, y - h - 0.010, 0.350, branch_y + 0.010)
  draw_arrow(xs[5] + w / 2, y - h - 0.010, 0.672, branch_y + 0.010)

  grid.roundrect(
    x = unit(0.500, "npc"),
    y = unit(0.095, "npc"),
    width = unit(0.865, "npc"),
    height = unit(0.080, "npc"),
    r = unit(0.016, "npc"),
    gp = gpar(fill = palette$pale, col = palette$line, lwd = 0.9)
  )
  text_grob(
    wrap_text("Reading rule: do not compare these as one shrinking group. Each number answers a different question: clinical visit inclusion, assembly QC, selected profiles, consecutive visit pairs, or Not_UTI -> UTI transitions.", 145),
    0.085, 0.120, 0.80, size = 8.6, color = palette$muted
  )
  text_grob("Sources: table_01_cohort_episode_flow.csv; pipeline_denominator_summary.csv; vf_longitudinal_transitions.csv; vf_transition_case_index.csv",
            0.045, 0.035, 0.9, size = 7.4, color = palette$faint)
}

draw_funnel_band <- function(y_top, h, top_w, bottom_w, fill, line,
                             title, metric, detail, idx, accent) {
  x_center <- 0.410
  x_top_l <- x_center - top_w / 2
  x_top_r <- x_center + top_w / 2
  x_bot_l <- x_center - bottom_w / 2
  x_bot_r <- x_center + bottom_w / 2
  y_bottom <- y_top - h
  grid.polygon(
    x = unit(c(x_top_l, x_top_r, x_bot_r, x_bot_l), "npc"),
    y = unit(c(y_top, y_top, y_bottom, y_bottom), "npc"),
    gp = gpar(fill = fill, col = line, lwd = 1.2)
  )
  grid.roundrect(
    x = unit(x_top_l + 0.037, "npc"),
    y = unit(y_top - 0.038, "npc"),
    width = unit(0.036, "npc"),
    height = unit(0.046, "npc"),
    r = unit(0.006, "npc"),
    gp = gpar(fill = accent, col = NA)
  )
  grid.text(idx, x = unit(x_top_l + 0.037, "npc"), y = unit(y_top - 0.038, "npc"),
            gp = gpar(fontsize = 8.2, col = palette$white, fontface = "bold"))
  text_grob(title, x_top_l + 0.070, y_top - 0.022, top_w - 0.100, size = 10.5, bold = TRUE)
  text_grob(metric, x_center, y_top - 0.077, top_w, size = 22, color = accent,
            bold = TRUE, just = c("center", "top"))
  text_grob(wrap_text(detail, 58), x_center, y_top - 0.130, top_w * 0.82, size = 8.2,
            color = palette$muted, just = c("center", "top"))
}

draw_breakdown_table <- function() {
  x <- 0.705
  y <- 0.735
  w <- 0.245
  h <- 0.395
  grid.roundrect(
    x = unit(x + w / 2, "npc"),
    y = unit(y - h / 2, "npc"),
    width = unit(w, "npc"),
    height = unit(h, "npc"),
    r = unit(0.016, "npc"),
    gp = gpar(fill = palette$pale, col = palette$line, lwd = 1)
  )
  text_grob("All status comparisons", x + 0.022, y - 0.030, w - 0.044, size = 11, bold = TRUE)
  text_grob("between consecutive clinical visits", x + 0.022, y - 0.064, w - 0.044, size = 7.6, color = palette$faint)

  order <- c("Not_UTI->Not_UTI", "UTI->Not_UTI", "Not_UTI->UTI", "UTI->UTI")
  fills <- c(palette$blue_fill, palette$slate_fill, palette$orange_fill, palette$yellow_fill)
  accents <- c(palette$blue, palette$faint, palette$orange, palette$yellow)
  for (i in seq_along(order)) {
    yy <- y - 0.112 - (i - 1) * 0.066
    n <- as.integer(transition_breakdown[[order[i]]])
    grid.roundrect(
      x = unit(x + 0.030, "npc"),
      y = unit(yy, "npc"),
      width = unit(0.020, "npc"),
      height = unit(0.026, "npc"),
      r = unit(0.003, "npc"),
      gp = gpar(fill = fills[i], col = accents[i], lwd = 0.7)
    )
    text_grob(order[i], x + 0.055, yy + 0.014, 0.125, size = 8.2, color = palette$muted)
    text_grob(fmt_n(n), x + w - 0.030, yy + 0.014, 0.035, size = 9.0, color = accents[i],
              bold = TRUE, just = c("right", "top"))
  }
  text_grob(
    wrap_text(
      paste0("Only the ", fmt_n(counts$clinical_not_uti_to_uti),
             " Not_UTI -> UTI transitions enter the focused review."),
      39
    ),
    x + 0.022, y - 0.342, w - 0.044, size = 7.8, color = palette$muted
  )
}

draw_funnel <- function() {
  draw_title(
    "Transition funnel: status comparisons to linked Not_UTI -> UTI transitions",
    "This companion view separates the transition question from the broader clinical-to-genomics workflow."
  )

  draw_funnel_band(
    0.790, 0.178, 0.610, 0.500, palette$blue_fill, "#BAE6FD",
    "All consecutive status comparisons",
    paste0(fmt_n(counts$total_transitions), " comparisons"),
    paste0(fmt_n(counts$transition_participants), " participants with consecutive clinical visits"),
    "1", palette$blue
  )
  draw_funnel_band(
    0.552, 0.165, 0.430, 0.350, palette$orange_fill, "#FED7AA",
    "Not_UTI -> UTI transitions",
    paste0(fmt_n(counts$clinical_not_uti_to_uti), " transitions"),
    "Status comparisons where a Not_UTI visit is followed by a UTI visit",
    "2", palette$orange
  )
  draw_funnel_band(
    0.335, 0.150, 0.310, 0.250, palette$green_fill, "#BBF7D0",
    "WGS/VF-linked transitions",
    paste0(fmt_n(counts$linked_not_uti_to_uti), " transitions"),
    paste0(fmt_n(counts$missing_linked_not_uti_to_uti),
           " transition is missing a usable VF-ready endpoint"),
    "3", palette$green
  )

  draw_arrow(0.410, 0.595, 0.410, 0.565)
  draw_arrow(0.410, 0.372, 0.410, 0.348)
  draw_breakdown_table()

  grid.roundrect(
    x = unit(0.500, "npc"),
    y = unit(0.085, "npc"),
    width = unit(0.865, "npc"),
    height = unit(0.080, "npc"),
    r = unit(0.016, "npc"),
    gp = gpar(fill = palette$pale, col = palette$line, lwd = 0.9)
  )
  text_grob(
    "Interpretation boundary: these transitions are evidence buckets for follow-up, not proven mechanisms. The linked subset is smaller because both endpoint visits must have usable WGS/VF evidence.",
    0.085, 0.110, 0.80, size = 8.6, color = palette$muted
  )
  text_grob("Source: results/vf/vf_transition_case_index.csv",
            0.045, 0.035, 0.9, size = 7.4, color = palette$faint)
}

svg_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

svg_width <- 2600
svg_height <- 1463
sx <- function(x) round(x * svg_width, 1)
sy <- function(y) round((1 - y) * svg_height, 1)
sw <- function(w) round(w * svg_width, 1)
sh <- function(h) round(h * svg_height, 1)

svg_text <- function(label, x, y, size = 24, fill = palette$ink,
                     weight = "400", anchor = "start",
                     wrap = NULL, line_height = 1.12) {
  lines <- if (is.null(wrap)) label else strwrap(label, width = wrap)
  if (length(lines) == 0) lines <- ""
  x_px <- sx(x)
  y_px <- sy(y) + size
  tspans <- character()
  for (i in seq_along(lines)) {
    dy <- if (i == 1) 0 else round(size * line_height, 1)
    tspans <- c(
      tspans,
      sprintf(
        '<tspan x="%.1f" dy="%.1f">%s</tspan>',
        x_px, dy, svg_escape(lines[i])
      )
    )
  }
  sprintf(
    '<text x="%.1f" y="%.1f" font-family="Arial, Helvetica, sans-serif" font-size="%.1f" fill="%s" font-weight="%s" text-anchor="%s">%s</text>',
    x_px, y_px, size, fill, weight, anchor, paste(tspans, collapse = "")
  )
}

svg_roundrect <- function(x, y, w, h, fill, stroke = palette$line,
                          stroke_width = 2, radius = 30) {
  sprintf(
    '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%.1f" fill="%s" stroke="%s" stroke-width="%.1f"/>',
    sx(x), sy(y), sw(w), sh(h), radius, fill, stroke, stroke_width
  )
}

svg_line <- function(x1, y1, x2, y2, color = "#94A3B8",
                     width = 4, marker = TRUE) {
  marker_attr <- if (marker) ' marker-end="url(#arrowhead)"' else ""
  sprintf(
    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.1f" stroke-linecap="round"%s/>',
    sx(x1), sy(y1), sx(x2), sy(y2), color, width, marker_attr
  )
}

svg_card <- function(x, y, w, h, idx, title, metric, detail,
                     fill = palette$slate_fill, line = palette$line,
                     accent = palette$ink) {
  c(
    svg_roundrect(x, y, w, h, fill, line, 2.4, 32),
    svg_roundrect(x + 0.006, y - 0.014, 0.036, 0.046, accent, accent, 0, 8),
    svg_text(idx, x + 0.024, y - 0.049, 17, palette$white, "700", "middle"),
    svg_text(title, x + 0.050, y - 0.038, 21, palette$ink, "700", wrap = 23, line_height = 0.95),
    svg_text(metric, x + 0.020, y - 0.118, 38, accent, "700"),
    svg_text(detail, x + 0.020, y - 0.170, 17, palette$faint, "400", wrap = 36, line_height = 1.08)
  )
}

svg_note <- function(x, y, label) {
  svg_text(label, x, y, 16.5, palette$faint, "400", wrap = 48, line_height = 1.08)
}

svg_title <- function(title, subtitle) {
  c(
    sprintf('<rect x="0" y="0" width="%d" height="%d" fill="%s"/>',
            svg_width, svg_height, palette$white),
    svg_text(title, 0.045, 0.955, 48, palette$ink, "700"),
    svg_text(subtitle, 0.045, 0.905, 23, palette$faint, "400", wrap = 135),
    svg_line(0.045, 0.855, 0.955, 0.855, palette$line, 2.2, marker = FALSE)
  )
}

svg_defs <- c(
  '<defs>',
  '  <marker id="arrowhead" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">',
  '    <path d="M 0 0 L 10 4 L 0 8 z" fill="#94A3B8"/>',
  '  </marker>',
  '</defs>'
)

write_svg <- function(path, body) {
  lines <- c(
    sprintf(
      '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" role="img">',
      svg_width, svg_height, svg_width, svg_height
    ),
    svg_defs,
    body,
    '</svg>'
  )
  writeLines(lines, path, useBytes = TRUE)
}

write_ladder_svg <- function(path) {
  w <- 0.168
  h <- 0.194
  y <- 0.760
  xs <- c(0.040, 0.232, 0.424, 0.616, 0.808)
  centers_y <- y - h / 2

  body <- c(
    svg_title(
      "Numbered workflow denominators with plain-language units",
      "Each card title names the unit; the large number shows the count for that stage."
    ),
    svg_card(xs[1], y, w, h, "1", "Classified clinical visits",
             paste0(fmt_n(counts$classified_episodes), " visits"),
             paste0(fmt_n(counts$classified_participants), " participants; ",
                    fmt_n(counts$classified_uti), " UTI, ",
                    fmt_n(counts$classified_not_uti), " Not_UTI before exclusions"),
             fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange),
    svg_card(xs[2], y, w, h, "2", "Included clinical visits",
             paste0(fmt_n(counts$included_episodes), " visits"),
             paste0(fmt_n(counts$included_participants), " participants; ",
                    fmt_n(counts$included_uti), " UTI, ",
                    fmt_n(counts$included_not_uti), " Not_UTI; 2 exclusions"),
             fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange),
    svg_card(xs[3], y, w, h, "3", "Assembly candidates for QC",
             fmt_n(counts$assembly_records),
             "Assembler alternatives stay explicit here; these are not extra visits",
             fill = palette$blue_fill, line = "#BAE6FD", accent = palette$blue),
    svg_card(xs[4], y, w, h, "4", "Selected genomic profiles",
             paste0(fmt_n(counts$canonical_episodes), " profiles"),
             paste0(fmt_n(counts$canonical_participants),
                    " participants; VF calls; one profile per usable visit"),
             fill = palette$blue_fill, line = "#BAE6FD", accent = palette$blue),
    svg_card(xs[5], y, w, h, "5", "VF/model-ready profiles",
             paste0(fmt_n(counts$vf_ready_episodes), " profiles"),
             paste0(fmt_n(counts$vf_ready_uti), " UTI, ",
                    fmt_n(counts$vf_ready_not_uti),
                    " Not_UTI; 27 included visits lack VF-ready evidence"),
             fill = palette$green_fill, line = "#BBF7D0", accent = palette$green)
  )
  for (i in 1:4) {
    body <- c(body, svg_line(xs[i] + w + 0.003, centers_y, xs[i + 1] - 0.003, centers_y))
  }
  body <- c(
    body,
    svg_card(0.222, 0.335, 0.256, 0.174, "6", "Longitudinal visit subset",
             paste0(fmt_n(counts$longitudinal_episodes), " visits"),
             paste0(fmt_n(counts$longitudinal_participants), " participants; ",
                    fmt_n(counts$longitudinal_comparisons), " consecutive visit pairs"),
             fill = palette$slate_fill, line = palette$line, accent = palette$blue),
    svg_card(0.544, 0.335, 0.256, 0.174, "7", "Not_UTI -> UTI transitions",
             paste0(fmt_n(counts$clinical_not_uti_to_uti), " transitions"),
             paste0(fmt_n(counts$linked_not_uti_to_uti),
                    " WGS/VF-linked; ",
                    fmt_n(counts$missing_linked_not_uti_to_uti), " missing endpoint"),
             fill = palette$slate_fill, line = palette$line, accent = palette$orange),
    svg_line(xs[5] + w / 2, y - h - 0.010, 0.350, 0.335 + 0.010),
    svg_line(xs[5] + w / 2, y - h - 0.010, 0.672, 0.335 + 0.010),
    svg_roundrect(0.0675, 0.135, 0.865, 0.080, palette$pale, palette$line, 1.8, 24),
    svg_text("Reading rule: do not compare these as one shrinking group. Each number answers a different question: clinical visit inclusion, assembly QC, selected profiles, consecutive visit pairs, or Not_UTI -> UTI transitions.",
             0.085, 0.120, 19.2, palette$muted, "400", wrap = 122),
    svg_text("Sources: table_01_cohort_episode_flow.csv; pipeline_denominator_summary.csv; vf_longitudinal_transitions.csv; vf_transition_case_index.csv",
             0.045, 0.035, 16.4, palette$faint, "400")
  )
  write_svg(path, body)
}

svg_label_box <- function(label, x, y, w, h = 0.040) {
  c(
    svg_roundrect(x, y, w, h, palette$white, palette$line, 1.3, 12),
    svg_text(label, x + 0.010, y - 0.020, 15.5, palette$muted, "400", wrap = 45)
  )
}

write_unit_story_svg <- function(path) {
  body <- c(
    svg_title(
      "Why the counts change: each stage counts a different thing",
      "Clinical visits expand into assembly candidates for QC, then collapse back to one selected genomic profile per usable visit."
    ),
    svg_text("1. Clinical curation", 0.055, 0.729, 23, palette$orange, "700"),
    svg_card(0.230, 0.770, 0.220, 0.135, "1", "Classified clinical visits",
             paste0(fmt_n(counts$classified_episodes), " visits"),
             "participant sample visits",
             fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange),
    svg_card(0.555, 0.770, 0.240, 0.135, "2", "Included clinical visits",
             paste0(fmt_n(counts$included_episodes), " visits"),
             paste0(fmt_n(counts$included_uti), " UTI, ",
                    fmt_n(counts$included_not_uti), " Not_UTI"),
             fill = palette$orange_fill, line = "#FED7AA", accent = palette$orange),
    svg_line(0.455, 0.705, 0.548, 0.705),
    svg_label_box("-2 manual exclusions", 0.455, 0.755, 0.110),
    svg_text("2. Temporary assembly inventory", 0.055, 0.539, 23, palette$blue, "700"),
    svg_card(0.370, 0.575, 0.300, 0.135, "3", "Assembly candidates / QC",
             paste0(fmt_n(counts$assembly_records), " candidates"),
             "genome files, not extra visits",
             fill = palette$blue_fill, line = "#BAE6FD", accent = palette$blue),
    svg_line(0.675, 0.635, 0.675, 0.585),
    svg_label_box("counting changes: visits -> assembly candidates", 0.695, 0.642, 0.220),
    svg_text("3. Return to one profile per visit", 0.055, 0.369, 23, palette$green, "700"),
    svg_card(0.370, 0.405, 0.300, 0.135, "4", "Selected genomic profiles",
             paste0(fmt_n(counts$canonical_episodes), " profiles"),
             "one profile per usable visit",
             fill = palette$green_fill, line = "#BBF7D0", accent = palette$green),
    svg_line(0.520, 0.440, 0.520, 0.405),
    svg_label_box("QC selects one profile per usable visit", 0.695, 0.455, 0.230),
    svg_text("4. Analysis-specific layers", 0.055, 0.207, 23, palette$muted, "700"),
    svg_card(0.215, 0.250, 0.220, 0.115, "5A", "VF/model-ready table",
             paste0(fmt_n(counts$vf_ready_episodes), " profiles"),
             paste0(fmt_n(counts$vf_ready_uti), " UTI, ",
                    fmt_n(counts$vf_ready_not_uti), " Not_UTI"),
             fill = palette$slate_fill, line = palette$line, accent = palette$green),
    svg_card(0.475, 0.250, 0.220, 0.115, "5B", "Longitudinal visit layer",
             paste0(fmt_n(counts$longitudinal_episodes), " visits"),
             paste0(fmt_n(counts$longitudinal_comparisons), " consecutive visit pairs"),
             fill = palette$slate_fill, line = palette$line, accent = palette$blue),
    svg_card(0.735, 0.250, 0.220, 0.115, "5C", "Focused transition layer",
             paste0(fmt_n(counts$clinical_not_uti_to_uti), " transitions"),
             paste0(fmt_n(counts$linked_not_uti_to_uti), " WGS/VF-linked"),
             fill = palette$slate_fill, line = palette$line, accent = palette$orange),
    svg_line(0.520, 0.270, 0.325, 0.252),
    svg_line(0.520, 0.270, 0.585, 0.252),
    svg_line(0.520, 0.270, 0.845, 0.252),
    svg_roundrect(0.0675, 0.095, 0.865, 0.060, palette$pale, palette$line, 1.8, 24),
    svg_text("Mental model: 585/583 are clinical visits; 1,291 is a temporary list of assembly candidates; 556 is selected genomic profiles; 538, 394, 11, and 10 are downstream analysis-specific subsets.",
             0.085, 0.081, 17.5, palette$muted, "400", wrap = 132),
    svg_text("Sources: pipeline_denominator_summary.csv; table_01_cohort_episode_flow.csv; vf_longitudinal_transitions.csv; vf_transition_case_index.csv",
             0.045, 0.018, 15.5, palette$faint, "400")
  )
  write_svg(path, body)
}

svg_poly <- function(points, fill, stroke, stroke_width = 2.4) {
  pts <- vapply(
    points,
    function(p) paste0(sx(p[1]), ",", sy(p[2])),
    character(1)
  )
  sprintf(
    '<polygon points="%s" fill="%s" stroke="%s" stroke-width="%.1f"/>',
    paste(pts, collapse = " "), fill, stroke, stroke_width
  )
}

svg_funnel_band <- function(y_top, h, top_w, bottom_w, fill, line,
                            title, metric, detail, idx, accent) {
  x_center <- 0.410
  x_top_l <- x_center - top_w / 2
  x_top_r <- x_center + top_w / 2
  x_bot_l <- x_center - bottom_w / 2
  x_bot_r <- x_center + bottom_w / 2
  y_bottom <- y_top - h
  c(
    svg_poly(
      list(c(x_top_l, y_top), c(x_top_r, y_top), c(x_bot_r, y_bottom), c(x_bot_l, y_bottom)),
      fill, line
    ),
    svg_roundrect(x_top_l + 0.019, y_top - 0.014, 0.036, 0.046, accent, accent, 0, 8),
    svg_text(idx, x_top_l + 0.037, y_top - 0.049, 17, palette$white, "700", "middle"),
    svg_text(title, x_top_l + 0.070, y_top - 0.044, 24, palette$ink, "700"),
    svg_text(metric, x_center, y_top - 0.118, 50, accent, "700", "middle"),
    svg_text(detail, x_center, y_top - 0.174, 18.5, palette$muted, "400", "middle", wrap = 62)
  )
}

svg_breakdown_table <- function() {
  x <- 0.705
  y <- 0.735
  w <- 0.245
  h <- 0.395
  order <- c("Not_UTI->Not_UTI", "UTI->Not_UTI", "Not_UTI->UTI", "UTI->UTI")
  fills <- c(palette$blue_fill, palette$slate_fill, palette$orange_fill, palette$yellow_fill)
  accents <- c(palette$blue, palette$faint, palette$orange, palette$yellow)

  body <- c(
    svg_roundrect(x, y, w, h, palette$pale, palette$line, 2, 24),
    svg_text("All status comparisons", x + 0.022, y - 0.052, 24, palette$ink, "700"),
    svg_text("between consecutive clinical visits", x + 0.022, y - 0.091, 17, palette$faint, "400")
  )
  for (i in seq_along(order)) {
    yy <- y - 0.112 - (i - 1) * 0.066
    n <- as.integer(transition_breakdown[[order[i]]])
    body <- c(
      body,
      svg_roundrect(x + 0.020, yy + 0.013, 0.020, 0.026, fills[i], accents[i], 1.2, 5),
      svg_text(order[i], x + 0.055, yy - 0.006, 18.5, palette$muted, "400"),
      svg_text(fmt_n(n), x + w - 0.030, yy - 0.006, 20, accents[i], "700", "end")
    )
  }
  c(
    body,
    svg_text(
      paste0("Only the ", fmt_n(counts$clinical_not_uti_to_uti),
             " Not_UTI -> UTI transitions enter the focused review."),
      x + 0.022, y - 0.368, 17.2, palette$muted, "400", wrap = 34
    )
  )
}

write_funnel_svg <- function(path) {
  body <- c(
    svg_title(
      "Transition funnel: status comparisons to linked Not_UTI -> UTI transitions",
      "This companion view separates the transition question from the broader clinical-to-genomics workflow."
    ),
    svg_funnel_band(
      0.790, 0.178, 0.610, 0.500, palette$blue_fill, "#BAE6FD",
      "All consecutive status comparisons",
      paste0(fmt_n(counts$total_transitions), " comparisons"),
      paste0(fmt_n(counts$transition_participants), " participants with consecutive clinical visits"),
      "1", palette$blue
    ),
    svg_funnel_band(
      0.552, 0.165, 0.430, 0.350, palette$orange_fill, "#FED7AA",
      "Not_UTI -> UTI transitions",
      paste0(fmt_n(counts$clinical_not_uti_to_uti), " transitions"),
      "Status comparisons where a Not_UTI visit is followed by a UTI visit",
      "2", palette$orange
    ),
    svg_funnel_band(
      0.335, 0.150, 0.310, 0.250, palette$green_fill, "#BBF7D0",
      "WGS/VF-linked transitions",
      paste0(fmt_n(counts$linked_not_uti_to_uti), " transitions"),
      paste0(fmt_n(counts$missing_linked_not_uti_to_uti),
             " transition is missing a usable VF-ready endpoint"),
      "3", palette$green
    ),
    svg_line(0.410, 0.595, 0.410, 0.565),
    svg_line(0.410, 0.372, 0.410, 0.348),
    svg_breakdown_table(),
    svg_roundrect(0.0675, 0.125, 0.865, 0.080, palette$pale, palette$line, 1.8, 24),
    svg_text("Interpretation boundary: these transitions are evidence buckets for follow-up, not proven mechanisms. The linked subset is smaller because both endpoint visits must have usable WGS/VF evidence.",
             0.085, 0.110, 19.2, palette$muted, "400", wrap = 120),
    svg_text("Source: results/vf/vf_transition_case_index.csv",
             0.045, 0.035, 16.4, palette$faint, "400")
  )
  write_svg(path, body)
}

render_figure <- function(path, draw_fun, type = c("png", "svg")) {
  type <- match.arg(type)
  if (type == "png") {
    png(path, width = 2600, height = 1463, res = 220, bg = "white")
  } else {
    svg(path, width = 12.0, height = 6.75, bg = "white")
  }
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  draw_fun()
}

story_png <- file.path(paths$out_dir, "08_unit_aware_denominator_story.png")
story_svg <- file.path(paths$out_dir, "08_unit_aware_denominator_story.svg")
ladder_png <- file.path(paths$out_dir, "09_case_count_workflow_ladder.png")
ladder_svg <- file.path(paths$out_dir, "09_case_count_workflow_ladder.svg")
funnel_png <- file.path(paths$out_dir, "10_transition_count_funnel.png")
funnel_svg <- file.path(paths$out_dir, "10_transition_count_funnel.svg")

render_figure(story_png, draw_unit_story, "png")
render_figure(ladder_png, draw_ladder, "png")
render_figure(funnel_png, draw_funnel, "png")
write_unit_story_svg(story_svg)
write_ladder_svg(ladder_svg)
write_funnel_svg(funnel_svg)

md_lines <- c(
  "# Workflow Denominator Map",
  "",
  "This guide explains why the numbers change as the project moves from clinical visits to genome assemblies, selected genomic profiles, and transition-focused comparisons. The point is that the workflow is not one shrinking cohort; it temporarily changes what is being counted.",
  "",
  "## Plain-Language Terms",
  "",
  "| Term | Meaning | Count examples |",
  "| :--- | :--- | :--- |",
  "| Clinical visit | One participant urine-sample visit/timepoint with a UTI or Not_UTI status. | 585 classified; 583 included |",
  "| Assembly candidate | One genome assembly file evaluated during QC. A single clinical visit can have more than one candidate assembly. | 1,291 |",
  "| Selected genomic profile | The one chosen genome assembly plus VF profile for a usable clinical visit. | 556 |",
  "| Consecutive VF visit pair | Two neighboring selected genomic profiles from the same participant. | 394 |",
  "| Not_UTI -> UTI transition | A consecutive visit comparison where a Not_UTI visit is followed by a UTI visit. | 11 total; 10 WGS/VF-linked |",
  "",
  "## What Is Going On?",
  "",
  "The confusing part is `1,291`. That number is not a larger set of people or clinical visits. It is a temporary list of assembly candidates where assembler alternatives and FASTA files are kept separately for QC. After QC, the workflow collapses back to one selected genomic profile per usable clinical visit, which is why the profile count becomes `556`.",
  "",
  "A cleaner way to read the main story is: `585` classified clinical visits -> `583` included clinical visits -> temporary expansion to `1,291` assembly candidates -> `556` selected genomic profiles -> analysis-specific subsets such as `538` longitudinal visits, `394` consecutive VF visit pairs, `11` Not_UTI -> UTI transitions, and `10` WGS/VF-linked transitions.",
  "",
  "## Images",
  "",
  "![Unit-aware denominator story](figures/workflow_flowchart/08_unit_aware_denominator_story.png)",
  "",
  "[SVG version](figures/workflow_flowchart/08_unit_aware_denominator_story.svg)",
  "",
  "![Numbered denominator ladder](figures/workflow_flowchart/09_case_count_workflow_ladder.png)",
  "",
  "[SVG version](figures/workflow_flowchart/09_case_count_workflow_ladder.svg)",
  "",
  "![Transition funnel](figures/workflow_flowchart/10_transition_count_funnel.png)",
  "",
  "[SVG version](figures/workflow_flowchart/10_transition_count_funnel.svg)",
  "",
  "## How To Read The Counts",
  "",
  "- `583` is the included clinical-visit denominator after manual curation exclusions.",
  "- `1,291` is larger because assembly QC counts genome assembly candidates, not extra people or visits.",
  "- `556` is the selected genomic profile denominator: one selected QC-pass assembly and VF profile per usable clinical visit.",
  "- `538` is the repeated-measures longitudinal visit subset; it produces `394` consecutive within-participant VF visit pairs.",
  "- `11` is the Not_UTI -> UTI transition count; `10` of those transitions have WGS/VF-linked endpoints.",
  "",
  "## Source Table",
  "",
  "| Step | Count shown | Human-readable unit | Technical source unit | Source | Notes |",
  "| :--- | :--- | :--- | :--- | :--- | :--- |",
  paste0("| Classified clinical visits | ", fmt_n(counts$classified_episodes), " | clinical visit | clinical_episode | `results/qc/pipeline_denominator_summary.csv` | Before primary manual curation exclusions; ", fmt_n(counts$classified_uti), " UTI and ", fmt_n(counts$classified_not_uti), " Not_UTI. |"),
  paste0("| Included clinical visits | ", fmt_n(counts$included_episodes), " | clinical visit | clinical_episode | `results/qc/pipeline_denominator_summary.csv` | ", fmt_n(counts$included_participants), " participants; ", fmt_n(counts$included_uti), " UTI and ", fmt_n(counts$included_not_uti), " Not_UTI. |"),
  paste0("| Assembly candidates for QC | ", fmt_n(counts$assembly_records), " | assembly candidate | assembly | `results/qc/pipeline_denominator_summary.csv` | Includes assembler alternatives, so it is not a clinical-visit count. |"),
  paste0("| Selected genomic profiles | ", fmt_n(counts$canonical_episodes), " | selected profile | participant_timepoint | `results/qc/pipeline_denominator_summary.csv` | ", fmt_n(counts$canonical_participants), " participants; one selected QC-pass profile per usable clinical visit. |"),
  paste0("| VF/model-ready profiles | ", fmt_n(counts$vf_ready_episodes), " | selected profile | participant_timepoint | `results/qc/pipeline_denominator_summary.csv` | ", fmt_n(counts$vf_ready_uti), " UTI and ", fmt_n(counts$vf_ready_not_uti), " Not_UTI. |"),
  paste0("| Longitudinal visit subset | ", fmt_n(counts$longitudinal_episodes), " | clinical visit with selected profile | participant_timepoint | `results/summary/table_01_cohort_episode_flow.csv` | ", fmt_n(counts$longitudinal_participants), " participants represented in VF transition output. |"),
  paste0("| Consecutive VF visit pairs | ", fmt_n(counts$longitudinal_comparisons), " | consecutive visit pair | participant_timepoint pair | `results/vf/vf_longitudinal_transitions.csv` | Derived from the longitudinal subset; ", fmt_n(counts$longitudinal_comparison_participants), " participants. |"),
  paste0("| All consecutive status comparisons | ", fmt_n(counts$total_transitions), " | status comparison | clinical_episode_transition | `results/vf/vf_transition_case_index.csv` | ", fmt_n(counts$transition_participants), " participants with consecutive ordered clinical states. |"),
  paste0("| Not_UTI -> UTI transitions | ", fmt_n(counts$clinical_not_uti_to_uti), " | transition | clinical_episode_transition | `results/vf/vf_transition_case_index.csv` | Focused transition denominator. |"),
  paste0("| WGS/VF-linked Not_UTI -> UTI transitions | ", fmt_n(counts$linked_not_uti_to_uti), " | linked transition | clinical_episode_transition | `results/vf/vf_transition_case_index.csv` | ", fmt_n(counts$missing_linked_not_uti_to_uti), " transition lacks a usable VF-ready endpoint. |"),
  "",
  "## Validation",
  "",
  "The figures and this Markdown file were generated by `scripts/create_workflow_case_count_flowchart.R`. The script re-reads the source CSVs and stops if the expected counts in this document drift from the current audited outputs."
)

writeLines(md_lines, paths$md, useBytes = TRUE)

message("Wrote:")
message("  ", paths$md)
message("  ", story_png)
message("  ", story_svg)
message("  ", ladder_png)
message("  ", ladder_svg)
message("  ", funnel_png)
message("  ", funnel_svg)

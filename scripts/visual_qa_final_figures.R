#!/usr/bin/env Rscript

# Render and transform every retained final figure for reproducible visual QA.
# Scientific figures are never edited in place: all QA derivatives are written
# below results/figure_audit/visual_qa/.

source("00_config.R")

suppressPackageStartupMessages({
  library(colorspace)
  library(dplyr)
  library(farver)
  library(png)
  library(readr)
  library(tidyr)
})

audit_dir <- file.path(DIR_RESULTS, "figure_audit")
manifest_path <- file.path(audit_dir, "final_figure_manifest.csv")
if (!file.exists(manifest_path)) stop("Missing final figure manifest: ", manifest_path, call. = FALSE)

manifest <- read_csv(manifest_path, show_col_types = FALSE)
if (nrow(manifest) != 18L || anyDuplicated(manifest$figure_id)) {
  stop("Visual QA requires the exact 18-family final figure manifest.", call. = FALSE)
}

qa_root <- file.path(audit_dir, "visual_qa")
qa_dirs <- c(
  pdf_render = file.path(qa_root, "pdf_render_150dpi"),
  a4 = file.path(qa_root, "a4_preview_1200px"),
  grey = file.path(qa_root, "greyscale"),
  deutan = file.path(qa_root, "deutan"),
  protan = file.path(qa_root, "protan")
)
invisible(lapply(qa_dirs, ensure_dir))

pdftoppm <- Sys.which("pdftoppm")
sips <- Sys.which("sips")
if (!nzchar(pdftoppm)) stop("pdftoppm is unavailable; PDF visual QA cannot proceed.", call. = FALSE)
if (!nzchar(sips)) stop("sips is unavailable; A4-size raster preview generation cannot proceed.", call. = FALSE)

srgb_to_linear <- function(x) {
  ifelse(x <= .03928, x / 12.92, ((x + .055) / 1.055)^2.4)
}

linear_to_srgb <- function(x) {
  x <- pmin(pmax(x, 0), 1)
  ifelse(x <= .03928 / 12.92, 12.92 * x, 1.055 * x^(1 / 2.4) - .055)
}

cvd_transform <- function(img, mode = c("deutan", "protan")) {
  mode <- match.arg(mode)
  h <- dim(img)[1]
  w <- dim(img)[2]
  rgb_matrix <- rbind(as.vector(img[, , 1]), as.vector(img[, , 2]), as.vector(img[, , 3]))
  transform <- if (mode == "deutan") {
    colorspace:::interpolate_cvd_transform(colorspace:::deutanomaly_cvd, 1)
  } else {
    colorspace:::interpolate_cvd_transform(colorspace:::protanomaly_cvd, 1)
  }
  mapped <- linear_to_srgb(transform %*% srgb_to_linear(rgb_matrix))
  out <- array(1, dim = c(h, w, if (dim(img)[3] == 4L) 4L else 3L))
  out[, , 1] <- matrix(mapped[1, ], nrow = h, ncol = w)
  out[, , 2] <- matrix(mapped[2, ], nrow = h, ncol = w)
  out[, , 3] <- matrix(mapped[3, ], nrow = h, ncol = w)
  if (dim(out)[3] == 4L) out[, , 4] <- img[, , 4]
  out
}

simulate_raster <- function(input, output, mode = c("grey", "deutan", "protan")) {
  mode <- match.arg(mode)
  img <- readPNG(input)
  if (length(dim(img)) != 3L || dim(img)[3] < 3L) stop("Unexpected PNG colour model: ", input, call. = FALSE)
  if (mode == "grey") {
    lum <- .2126 * img[, , 1] + .7152 * img[, , 2] + .0722 * img[, , 3]
    out <- array(1, dim = c(dim(img)[1], dim(img)[2], if (dim(img)[3] == 4L) 4L else 3L))
    out[, , 1] <- lum
    out[, , 2] <- lum
    out[, , 3] <- lum
    if (dim(out)[3] == 4L) out[, , 4] <- img[, , 4]
  } else {
    out <- cvd_transform(img, mode)
  }
  writePNG(out, output)
  invisible(output)
}

read_png_dimensions <- function(path) {
  x <- readPNG(path, info = TRUE)
  c(width_px = dim(x)[2], height_px = dim(x)[1])
}

qa_rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  id <- manifest$figure_id[[i]]
  png_path <- manifest$png_path[[i]]
  pdf_path <- manifest$pdf_path[[i]]
  if (!file.exists(png_path) || !file.exists(pdf_path)) stop("Missing final figure file for ", id, call. = FALSE)

  pdf_stub <- file.path(qa_dirs[["pdf_render"]], id)
  pdf_png <- paste0(pdf_stub, ".png")
  pdf_status <- system2(pdftoppm, c("-png", "-singlefile", "-r", "150", pdf_path, pdf_stub), stdout = FALSE, stderr = FALSE)
  if (!identical(pdf_status, 0L) || !file.exists(pdf_png) || file.info(pdf_png)$size <= 0) {
    stop("Poppler failed to render ", pdf_path, call. = FALSE)
  }

  a4_path <- file.path(qa_dirs[["a4"]], paste0(id, ".png"))
  a4_status <- system2(sips, c("-Z", "1200", png_path, "--out", a4_path), stdout = FALSE, stderr = FALSE)
  if (!identical(a4_status, 0L) || !file.exists(a4_path) || file.info(a4_path)$size <= 0) {
    stop("A4 preview generation failed for ", png_path, call. = FALSE)
  }

  grey_path <- file.path(qa_dirs[["grey"]], paste0(id, ".png"))
  deutan_path <- file.path(qa_dirs[["deutan"]], paste0(id, ".png"))
  protan_path <- file.path(qa_dirs[["protan"]], paste0(id, ".png"))
  simulate_raster(a4_path, grey_path, "grey")
  simulate_raster(a4_path, deutan_path, "deutan")
  simulate_raster(a4_path, protan_path, "protan")

  dims <- read_png_dimensions(a4_path)
  qa_rows[[i]] <- tibble(
    figure_id = id,
    figure_number = manifest$figure_number[[i]],
    class = manifest$figure_class[[i]],
    pdf_render_path = pdf_png,
    a4_preview_path = a4_path,
    greyscale_path = grey_path,
    deutan_path = deutan_path,
    protan_path = protan_path,
    preview_width_px = unname(dims[["width_px"]]),
    preview_height_px = unname(dims[["height_px"]]),
    all_derivatives_nonempty = all(file.info(c(pdf_png, a4_path, grey_path, deutan_path, protan_path))$size > 0)
  )
}
qa <- bind_rows(qa_rows)
if (nrow(qa) != 18L || any(!qa$all_derivatives_nonempty)) stop("Visual-QA derivative checks failed.", call. = FALSE)
write_csv(qa, file.path(audit_dir, "visual_qa_results.csv"))

palette <- tibble(
  semantic_level = c("UTI", "Not UTI", "Legacy Negative", "Unknown"),
  hex = c("#D55E00", "#0072B2", "#909090", "#CCCCCC")
) %>%
  mutate(
    contrast_on_white = as.numeric(colorspace::contrast_ratio(.data$hex, "#FFFFFF")),
    deutan_hex = colorspace::deutan(.data$hex),
    protan_hex = colorspace::protan(.data$hex)
  )
write_csv(palette, file.path(audit_dir, "palette_accessibility_checks.csv"))

pair_distance <- function(colours) {
  lab <- farver::decode_colour(colours, to = "lab")
  farver::compare_colour(lab[1, , drop = FALSE], lab[2, , drop = FALSE], from_space = "lab", method = "CIE2000")[[1]]
}
status_pair <- tibble(
  simulation = c("Original", "Deutan", "Protan"),
  uti = c(palette$hex[1], palette$deutan_hex[1], palette$protan_hex[1]),
  not_uti = c(palette$hex[2], palette$deutan_hex[2], palette$protan_hex[2])
) %>%
  rowwise() %>%
  mutate(ciede2000_distance = pair_distance(c(.data$uti, .data$not_uti))) %>%
  ungroup()
write_csv(status_pair, file.path(audit_dir, "operational_status_colour_separation.csv"))

writeLines(c(
  sprintf("Generated visual-QA derivatives for %d final figure families.", nrow(qa)),
  sprintf("Poppler PDF renders: %d/%d nonempty.", sum(file.info(qa$pdf_render_path)$size > 0), nrow(qa)),
  sprintf("A4-width previews: %d/%d nonempty.", sum(file.info(qa$a4_preview_path)$size > 0), nrow(qa)),
  sprintf("Greyscale simulations: %d/%d nonempty.", sum(file.info(qa$greyscale_path)$size > 0), nrow(qa)),
  sprintf("Deutan simulations: %d/%d nonempty.", sum(file.info(qa$deutan_path)$size > 0), nrow(qa)),
  sprintf("Protan simulations: %d/%d nonempty.", sum(file.info(qa$protan_path)$size > 0), nrow(qa)),
  "Automated derivative creation does not substitute for recorded human visual inspection."
), file.path(audit_dir, "visual_qa_generation_summary.txt"))

cat(sprintf("Visual-QA derivative generation passed: %d families, %d derivative files.\n", nrow(qa), 5L * nrow(qa)))

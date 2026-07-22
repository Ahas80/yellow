library(testthat)
source(here::here("R", "plot_helpers.R"))

test_that("canonical operational and legacy scales retain semantic colours", {
    expect_identical(
        uti_status_cols,
        c(UTI = "#D55E00", Not_UTI = "#0072B2", Unknown = "#CCCCCC")
    )
    expect_identical(
        infection_cols,
        c(UTI = "#D55E00", ASB = "#0072B2", Negative = "#909090",
          Unknown = "#CCCCCC")
    )
    expect_s3_class(scale_colour_uti_status(), "ScaleDiscrete")
    expect_s3_class(scale_fill_clinical_episode(), "ScaleDiscrete")
    display_scale <- scale_fill_operational_uti(reader_facing = TRUE)
    display_scale$train(c("UTI", "Not UTI"))
    expect_identical(unname(display_scale$map("Not UTI")), "#0072B2")
})

test_that("reader-facing recoders preserve missingness and reject bad levels", {
    status <- recode_operational_uti_status(
        c("UTI", "Not_UTI", "Unknown", NA), as_factor = FALSE
    )
    expect_identical(status, c("UTI", "Not UTI", "Unknown", NA_character_))
    expect_error(
        recode_operational_uti_status("ASB"),
        "unsupported level"
    )
    expect_identical(
        as.character(recode_clinical_episode_type(c("ASB", "Negative"))),
        c("ASB", "Negative")
    )
})

test_that("manual scale assertion catches palette mismatches", {
    expect_true(assert_ruti_scale_levels(
        factor(c("UTI", "Not_UTI"), levels = c("UTI", "Not_UTI")),
        uti_status_cols
    ))
    expect_error(
        assert_ruti_scale_levels(c("UTI", "ASB"), uti_status_cols),
        "ASB"
    )
    expect_error(
        assert_ruti_scale_levels(c("UTI", NA), uti_status_cols, allow_na = FALSE),
        "missing mapped values"
    )
})

test_that("deterministic helpers are stable and restore RNG state", {
    expect_identical(ruti_seed_from_key("Fig01"), ruti_seed_from_key("Fig01"))
    expect_false(identical(ruti_seed_from_key("Fig01"), ruti_seed_from_key("Fig02")))
    set.seed(42)
    before <- .Random.seed
    a <- with_ruti_seed(stats::runif(4), seed = 99)
    expect_identical(.Random.seed, before)
    b <- with_ruti_seed(stats::runif(4), seed = 99)
    expect_identical(a, b)
    expect_s3_class(position_jitter_ruti(key = "Fig01"), "PositionJitter")
})

test_that("stable case labels require a secret and do not expose identifiers", {
    ids <- c("participant-1", "participant-2", "participant-1", NA)
    salt <- "unit-test-secret-salt-2026"
    labels_a <- stable_case_labels(ids, salt = salt)
    labels_b <- stable_case_labels(rev(ids), salt = salt)
    expect_identical(labels_a[1], labels_a[3])
    expect_identical(labels_a[1], labels_b[4])
    expect_true(is.na(labels_a[4]))
    expect_false(any(grepl("participant", labels_a, fixed = TRUE), na.rm = TRUE))
    expect_error(stable_case_labels(ids, salt = "short"), "at least 16")
})

test_that("unsafe ratio estimates are classified and never plotted as finite", {
    prepared <- prepare_effect_estimates_for_plot(
        estimate = c(2, Inf, 0, 5, 1.5),
        conf_low = c(1.2, 1, 0, 0.5, 1),
        conf_high = c(3.4, Inf, 1, 20, 2),
        effect_scale = "ratio",
        separation = c(FALSE, FALSE, FALSE, TRUE, FALSE),
        limits = c(0.25, 8)
    )
    expect_identical(
        as.character(prepared$effect_status),
        c("estimable", "non_finite_estimate", "invalid_for_scale",
          "separation_or_instability", "estimable")
    )
    expect_true(all(is.na(prepared$plot_estimate[!prepared$effect_estimable])))
    expect_true(all(prepared$annotation_x[!prepared$effect_estimable] == 8))
    empty <- prepare_effect_estimates_for_plot(numeric(), effect_scale = "ratio")
    expect_identical(nrow(empty), 0L)
})

test_that("publication theme and save helper produce validated PNG/PDF metadata", {
    expect_s3_class(theme_ruti_publication(), "theme")
    out_dir <- tempfile("ruti_plot_test_")
    dir.create(out_dir)
    stem <- file.path(out_dir, "FigTest")
    manifest <- file.path(out_dir, "manifest.csv")
    p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                         ggplot2::aes(x, y)) +
        ggplot2::geom_point() +
        theme_ruti_publication()
    entry <- save_ruti_figure(
        p, stem, width = 3, height = 2, manifest_path = manifest,
        metadata = list(source_script = "tests/test_plot_helpers.R")
    )
    expect_true(file.exists(paste0(stem, ".png")))
    expect_true(file.exists(paste0(stem, ".pdf")))
    expect_true(file.exists(manifest))
    expect_gt(file.info(paste0(stem, ".png"))$size, 0)
    expect_gt(file.info(paste0(stem, ".pdf"))$size, 0)
    png_connection <- file(paste0(stem, ".png"), open = "rb")
    png_header <- readBin(png_connection, what = "raw", n = 24L)
    close(png_connection)
    png_width <- sum(as.integer(png_header[17:20]) * 256^(3:0))
    png_height <- sum(as.integer(png_header[21:24]) * 256^(3:0))
    expect_identical(c(png_height, png_width), c(600, 900))
    pdf_connection <- file(paste0(stem, ".pdf"), open = "rb")
    on.exit(close(pdf_connection), add = TRUE)
    expect_identical(readChar(pdf_connection, 5L, useBytes = TRUE), "%PDF-")
    expect_identical(entry$dpi, 300L)
    expect_identical(entry$background, "white")
    written <- utils::read.csv(manifest, stringsAsFactors = FALSE)
    expect_identical(written$figure_id, "FigTest")
    expect_identical(written$source_script, "tests/test_plot_helpers.R")
    expect_error(
        save_ruti_figure(p, file.path(out_dir, "bad"), 3, 2, dpi = 299,
                         write_manifest = FALSE),
        "at least 300"
    )
})

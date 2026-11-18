# Hi-C settings as W93 set them: the bin size never exceeds what plotgardener draws (bug 137), the
# colour scale's automatic top is the 95th percentile of the contacts in view (Juicebox's rule),
# rounded up (bug 203), the bin size drawn is shown under the colour key, and the settings' labels.

# ---- Bin size choice (no data needed) -------------------------------------------------------------

pick <- function(avail, width, ...) .pt_hic_pick_resolution(avail, width, list(...))

test_that("Auto never picks a bin size giving over 1,000 bins across the view (bug 137)", {
  avail <- c(1000, 5000, 10000, 25000)
  # 2 Mb: 1 kb bins would be 2,000 across the view (plotgardener would redraw at its own size);
  # 5 kb (400 bins) is the finest that fits, though it gives fewer than the 500 asked for.
  expect_identical(pick(avail, 2e6, hic_res_mode = "auto", hic_max_bins_per_view = 500), 5000)
  # Within the limit, unchanged: the coarsest size giving at least 500 bins.
  expect_identical(pick(avail, 5e6, hic_res_mode = "auto", hic_max_bins_per_view = 500), 10000)
  expect_identical(pick(avail, 1e6, hic_res_mode = "auto", hic_max_bins_per_view = 500), 1000)
  # A minimum above the limit can't be met: the finest size that fits.
  expect_identical(pick(avail, 5e6, hic_res_mode = "auto", hic_max_bins_per_view = 5000), 5000)
  # No minimum saved (an old settings file): 500.
  expect_identical(pick(avail, 5e6, hic_res_mode = "auto"), 10000)
})

test_that("Finest and Fixed are limited to 1,000 bins across the view too", {
  avail <- c(1000, 5000, 10000, 25000)
  expect_identical(pick(avail, 2e6, hic_res_mode = "finest"), 5000)
  expect_identical(pick(avail, 5e5, hic_res_mode = "finest"), 1000)
  expect_identical(pick(avail, 2e6, hic_res_mode = "fixed", hic_resolution = 1000), 5000)
  expect_identical(pick(avail, 2e6, hic_res_mode = "fixed", hic_resolution = 10000), 10000)
  expect_error(pick(avail, 2e6, hic_res_mode = "fixed", hic_resolution = 2000), "no 2,000 bp bins")
  # Even the coarsest size is too fine: the coarsest (plotgardener then picks its own).
  expect_identical(pick(c(1000, 5000), 2e7, hic_res_mode = "auto", hic_max_bins_per_view = 500), 5000)
})

# ---- Colour scale top and bin label -------------------------------------------------------------

test_that("the colour scale's top is rounded up to 4 significant digits, never below the value (bug 203)", {
  expect_gte(.pt_hic_round_top(0.00012345), 0.00012345)
  expect_equal(.pt_hic_round_top(0.00012345), 0.0001235)
  expect_equal(.pt_hic_round_top(450), 450)
  expect_equal(.pt_hic_round_top(2391.27), 2392)
  expect_equal(.pt_hic_round_top(1e-5), 1e-5)
  expect_equal(.pt_hic_round_top(6.57), 6.57)
  for (x in c(0.000999999, 123456.7, 0.3, 7e-9)) expect_gte(.pt_hic_round_top(x), x)
})

test_that("bin sizes are labelled as people write them", {
  expect_identical(.pt_hic_bin_label(500), "500 bp")
  expect_identical(.pt_hic_bin_label(5000), "5 kb")
  expect_identical(.pt_hic_bin_label(250000), "250 kb")
  expect_identical(.pt_hic_bin_label(2500000), "2.5 Mb")
})

test_that("with no Colour max, the colour scale tops at the 95th percentile of the contacts in view (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- list(chrom = "chr2", start = 100000000, end = 105000000, gene = NA_character_)
  settings <- test_settings("mouse")
  entry <- single_track("hic", test_rows("mouse")[["Ba"]], settings = settings)
  out <- settings_hic(list(), entry, settings, 7, region)
  # The contacts drawn (W95): Rectangle's strip.
  counts <- .pt_hic_drawn_counts(out$entry$data$matrix, "rectangle", region$start, region$end,
                                 out$entry$data$resolution, out$entry$height / 7)
  top <- out$dict[[.pt_shared_range_key("hic_zrange_", entry)]]
  expect_equal(top, unname(stats::quantile(counts, 0.95, na.rm = TRUE)))
  # Far below the maximum: the diagonal no longer takes the whole scale.
  expect_lt(top, max(counts, na.rm = TRUE) / 2)
})

test_that("the bin size drawn is shown under the colour key, and follows the Resolution setting (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  row <- test_rows("mouse")[["Ba"]]
  svg_of <- function(settings) {
    f <- tempfile(fileext = ".svg")
    on.exit(unlink(f), add = TRUE)
    suppressMessages(suppressWarnings(render_all(list(single_track("hic", row, settings = settings)), region, settings,
      device_fn = function(w, h) svglite::svglite(f, width = w, height = h))))
    grDevices::dev.off()
    paste(readLines(f, warn = FALSE), collapse = "\n")
  }
  auto <- svg_of(test_settings("mouse"))
  res <- load_hic(single_track("hic", row), region, test_settings("mouse"))$data$resolution
  expect_match(auto, paste0(">", .pt_hic_bin_label(res), " bins<"), fixed = TRUE)
  fixed <- svg_of(test_settings("mouse", list(hic_res_mode = "fixed", hic_resolution = 50000)))
  expect_match(fixed, ">50 kb bins<", fixed = TRUE)
})

test_that("Colour scale Log runs from the smallest value above 0 to the largest, and says so (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  row <- test_rows("mouse")[["Ba"]]
  lin_settings <- test_settings("mouse")
  log_settings <- test_settings("mouse", list(hic_colour_scale = "log"))
  out <- settings_hic(list(), single_track("hic", row, settings = log_settings), log_settings, 7, region)
  counts <- .pt_hic_drawn_counts(out$entry$data$matrix, "rectangle", region$start, region$end,
                                 out$entry$data$resolution, out$entry$height / 7)
  counts <- counts[is.finite(counts)]
  expect_true(out$entry$data$log_scale)
  expect_equal(out$dict[[.pt_shared_range_key("hic_zrange_", out$entry)]], max(counts))
  expect_equal(out$dict[[.pt_shared_range_key("hic_zlow_", out$entry)]], min(counts[counts > 0]))

  svg_of <- function(settings) {
    f <- tempfile(fileext = ".svg")
    on.exit(unlink(f), add = TRUE)
    suppressMessages(suppressWarnings(render_all(list(single_track("hic", row, settings = settings)), region, settings,
      device_fn = function(w, h) svglite::svglite(f, width = w, height = h))))
    grDevices::dev.off()
    paste(readLines(f, warn = FALSE), collapse = "\n")
  }
  lin <- svg_of(lin_settings)
  lg <- svg_of(log_settings)
  expect_false(grepl(">log scale<", lin, fixed = TRUE))
  expect_match(lg, ">log scale<", fixed = TRUE)
  fills <- function(svg) regmatches(svg, gregexpr("<polygon[^>]*fill: #[0-9A-Fa-f]+", svg))[[1]]
  expect_gt(length(fills(lg)), 0)
  # The transform itself, at one colour range: only colorTrans differs between the two.
  same <- list(zrange_cap = 20)
  expect_false(identical(fills(svg_of(modifyList(lin_settings, same))), fills(svg_of(modifyList(log_settings, same)))))
  # A Colour max with Log: the bottom stays the smallest value, the top is the max set.
  capped <- svg_of(test_settings("mouse", list(hic_colour_scale = "log", zrange_cap = 20)))
  expect_false(grepl("track failed", capped, fixed = TRUE))
  expect_match(capped, ">20<", fixed = TRUE)
})

test_that("a cleared Normalisation box keeps the saved normalisation instead of sending an empty one", {
  got <- .pt_collect_display_settings(list(s_hic_norm = ""), "s_", list(hic_norm = "VC"))
  expect_identical(got$hic_norm, "VC")
  expect_identical(.pt_collect_display_settings(list(s_hic_norm = " "), "s_", list())$hic_norm, "KR")
  expect_identical(.pt_collect_display_settings(list(s_hic_norm = "ICE"), "s_", list())$hic_norm, "ICE")
  # An empty one saved before (the old free-text box allowed it): KR.
  expect_identical(.pt_collect_display_settings(list(s_hic_norm = ""), "s_", list(hic_norm = ""))$hic_norm, "KR")
  expect_identical(.pt_render_cfg(list(hic_norm = " "), list(), tempdir())$hic_norm, "KR")
})

test_that("the colour key labels the colour range to 3 significant digits", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  row <- test_rows("mouse")[["Ba"]]
  settings <- test_settings("mouse", list(zrange_cap = 12.3456))
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  suppressMessages(suppressWarnings(render_all(list(single_track("hic", row, settings = settings)), region, settings,
    device_fn = function(w, h) svglite::svglite(f, width = w, height = h))))
  grDevices::dev.off()
  expect_match(paste(readLines(f, warn = FALSE), collapse = "\n"), ">12.3<", fixed = TRUE)
})

test_that("a Normalisation the file lacks says so and lists the ones it has (mouse and human datasets)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  settings <- test_settings("mouse", list(hic_norm = "BALANCE"))
  spec <- single_track("hic", test_rows("mouse")[["Ba"]], settings = settings)
  expect_error(load_hic(spec, region, settings), "^this file has no BALANCE normalisation at [0-9,]+ bp bins: in Figure settings, Hi-C, set Normalisation \\(hic_norm\\) to [A-Z_, ]+ or NONE[.]$")
  mcool <- file.path(DATA_ROOT, "hic_test", "LEUK_HEK_PJA27_inter_30.mcool")
  skip_if_not(file.exists(mcool))
  spec$path <- mcool
  wide <- list(chrom = "chr2", start = 100000000, end = 110000000, gene = NA_character_)
  msg <- tryCatch(load_hic(spec, wide, settings), error = conditionMessage)
  expect_identical(msg, "this file has no BALANCE normalisation: in Figure settings, Hi-C, set Normalisation (hic_norm) to KR, VC, VC_SQRT or NONE.")
  expect_lte(nchar(msg), 140)   # render_all() shows a message of up to 140 characters uncut, after "<track> track failed: "
})

# ---- Settings labels ------------------------------------------------------------------------------

test_that("the Hi-C settings say what they do: Min bins across view (Auto only), a normalisation list, Colour max", {
  html <- as.character(.pt_display_settings_ui(list(hic_norm = "weight"), "s_")$display)
  expect_false(grepl("Max bins per view", html, fixed = TRUE))
  expect_match(html, "Min bins across view:", fixed = TRUE)
  # Shown only when Resolution is Auto, the one mode it changes.
  expect_match(html, "data-display-if=\"input.s_hic_res_mode === &#39;auto&#39;\"", fixed = TRUE)
  expect_match(html, "Resolution:", fixed = TRUE)
  expect_match(html, "Colour max:", fixed = TRUE)
  expect_match(html, "Colour scale:", fixed = TRUE)
  expect_match(html, "<option value=\"log\"", fixed = TRUE)
  expect_match(html, "95th percentile", fixed = TRUE)
  # A list of the usual names, keeping a name the project already uses.
  for (n in c(HIC_NORMS, "weight")) expect_match(html, paste0("<option value=\"", n, "\""), fixed = TRUE)
  expect_match(html, "<option value=\"weight\" selected>", fixed = TRUE)
})

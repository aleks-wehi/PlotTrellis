# W95: the Hi-C Colour max is one setting, zrange_cap: NULL (Auto: the 95th percentile on a
# linear colour scale, the maximum on a log one), "max", "p95", "p90" or a number. The statistic
# is taken over the contacts the track draws (the Rectangle style's strip, Square's whole
# square), it round-trips through the settings panel, project_settings.json, samples.txt and the
# downloaded script, and the Hi-C normalisation errors say where to fix them.

# ---- The setting's values (no data needed) -----------------------------------------------------

test_that("a Colour max reads as Auto, a statistic or a number", {
  for (auto in list(NULL, NA, "", " ", "auto", "Auto", "NA", list())) expect_null(.pt_hic_zrange_cap(auto))
  expect_identical(.pt_hic_zrange_cap("max"), "max")
  expect_identical(.pt_hic_zrange_cap(" P95 "), "p95")
  expect_identical(.pt_hic_zrange_cap("p90"), "p90")
  expect_identical(.pt_hic_zrange_cap(50L), 50L)       # an old project's number
  expect_identical(.pt_hic_zrange_cap(12.5), 12.5)
  expect_identical(.pt_hic_zrange_cap("12.5"), 12.5)   # typed in the box, or from samples.txt
  for (bad in list("p99", "abc", Inf, c(0, 50), TRUE)) expect_true(is.na(.pt_hic_zrange_cap(bad)))
  expect_null(.pt_hic_zrange_cap("p99", bad = NULL))
})

# ---- The contacts drawn ------------------------------------------------------------------------

# A 1 kb view (1,000-2,000) at 100 bp bins, drawn 0.5 view widths high (Rectangle's default):
# the square's 55 contacts, valued 1 to 55, and contacts in the strip left and right of the
# square (value 0.5, and one 500).
SQUARE <- local({
  b <- seq(1000, 1900, by = 100)
  g <- expand.grid(x = b, y = b)
  g <- g[g$x <= g$y, ]
  data.frame(x = g$x, y = g$y, counts = seq_len(nrow(g)))
})
STRIP_ONLY <- data.frame(x = c(600, 700, 800, 900, 500, 1500, 1600),
                         y = c(1300, 1200, 1400, 1100, 1400, 2300, 2100),
                         counts = c(0.5, 0.5, 0.5, 0.5, 500, 0.5, 0.5))
# Read but not drawn: midpoint before the view, or further from the diagonal than the strip is high.
OUTSIDE <- data.frame(x = c(500, 1000, 2400), y = c(600, 2100, 2500), counts = c(1000, 1000, 1000))
MATRIX <- rbind(SQUARE, STRIP_ONLY, OUTSIDE)

test_that("Rectangle draws the strip around the view; Square both halves of its square; Triangle one", {
  rect <- .pt_hic_drawn_counts(MATRIX, "rectangle", 1000, 2000, 100, 0.5)
  expect_equal(sort(rect), sort(c(SQUARE$counts, STRIP_ONLY$counts)))
  sq <- .pt_hic_drawn_counts(SQUARE, "square", 1000, 2000, 100, 1)
  expect_equal(sort(sq), sort(c(SQUARE$counts, SQUARE$counts[SQUARE$x != SQUARE$y])))
  expect_equal(.pt_hic_drawn_counts(SQUARE, "triangle", 1000, 2000, 100, 0.5), SQUARE$counts)
  # A taller rectangle reaches further from the diagonal: the 1000 bp-apart contact at 2,100.
  expect_length(.pt_hic_drawn_counts(MATRIX, "rectangle", 1000, 2000, 100, 0.6), length(rect) + 1)
})

# settings_hic() with load_hic() replaced by the matrix above: the read it asks for and the top.
fake_hic <- function(style, cap = NULL, scale = "linear", height = NULL, matrix = MATRIX) {
  pads <- numeric(0)
  local_mocked_bindings(load_hic = function(entry, region_params, settings, pad = 0) {
    pads <<- c(pads, pad)
    entry$data <- list(matrix = if (style == "rectangle") matrix else matrix[matrix$x >= 1000 & matrix$y < 2000, ],
                       chrom = "chr1", resolution = 100)
    entry
  })
  entry <- list(viz = "hic", name = "HiC", hic_style = style, height = height)
  out <- settings_hic(list(), entry, list(zrange_cap = cap, hic_colour_scale = scale), 7,
                      list(chrom = "chr1", start = 1000, end = 2000))
  list(top = out$dict[["hic_zrange_HiC"]], low = out$dict[["hic_zlow_HiC"]], pad = pads)
}

test_that("each Colour max gives its statistic of the contacts drawn", {
  rect <- c(SQUARE$counts, STRIP_ONLY$counts)
  q <- function(x, p) unname(stats::quantile(x, p))
  expect_equal(fake_hic("rectangle", "max")$top, 500)
  expect_equal(fake_hic("rectangle", "p95")$top, q(rect, 0.95))
  expect_equal(fake_hic("rectangle", "p90")$top, q(rect, 0.90))
  expect_equal(fake_hic("rectangle", NULL)$top, q(rect, 0.95))           # Auto, Linear: 95th
  expect_equal(fake_hic("rectangle", NULL, "log")$top, 500)              # Auto, Log: max
  expect_equal(fake_hic("rectangle", NULL, "log")$low, 0.5)              # Log's bottom: drawn too
  expect_equal(fake_hic("rectangle", "p90", "log")$top, q(rect, 0.90))   # a choice holds on Log
  expect_identical(fake_hic("rectangle", 50)$top, list(override = c(0, 50)))
  expect_identical(fake_hic("rectangle", "50")$top, list(override = c(0, 50)))
  # The three statistics differ here, so a mix-up can't pass.
  expect_length(unique(c(500, q(rect, 0.95), q(rect, 0.90))), 3)
})

test_that("Rectangle's top comes from its strip, not the view's square", {
  sq_only <- SQUARE$counts
  rect <- fake_hic("rectangle", "p90")
  expect_equal(rect$top, unname(stats::quantile(c(sq_only, STRIP_ONLY$counts), 0.9)))
  expect_false(isTRUE(all.equal(rect$top, unname(stats::quantile(sq_only, 0.9)))))
  expect_false(isTRUE(all.equal(fake_hic("rectangle", "max")$top, max(sq_only))))
  # It reads the strip: as many bp beyond each end as the track is high, in view widths
  # (3.5 in high on a 7 in track: half the 1,000 bp view). Square and Triangle read the square.
  expect_equal(rect$pad, 500)
  expect_equal(fake_hic("rectangle", "p90", height = 7)$pad, 1000)
  expect_equal(fake_hic("square", "p90")$pad, 0)
  expect_equal(fake_hic("triangle", "p90")$pad, 0)
  # Square counts each contact off the diagonal twice (both halves drawn).
  expect_equal(fake_hic("square", "p90")$top,
               unname(stats::quantile(c(sq_only, SQUARE$counts[SQUARE$x != SQUARE$y]), 0.9)))
})

test_that("same-name Hi-C tracks share the largest top", {
  local_mocked_bindings(load_hic = function(entry, region_params, settings, pad = 0) {
    entry$data <- list(matrix = data.frame(x = 1000, y = 1000, counts = entry$value), chrom = "chr1",
                       resolution = 100)
    entry
  })
  dict <- list()
  for (v in c(3, 8, 5)) dict <- settings_hic(dict, list(viz = "hic", name = "HiC", value = v, height = 3.5),
                                            list(zrange_cap = "max"), 7,
                                            list(chrom = "chr1", start = 1000, end = 2000))$dict
  expect_equal(dict[["hic_zrange_HiC"]], 8)
  expect_identical(.pt_hic_colour_tops(c(dict, list(hic_zrange_other = list(override = c(0, 4))))),
                   c(hic_zrange_HiC = 8))
})

test_that("a Colour max that is not one says so in the track", {
  bad <- function(v) tryCatch(settings_hic(list(), list(viz = "hic", name = "HiC", height = 2), list(zrange_cap = v), 7,
                                           list(chrom = "chr1", start = 1, end = 2)), error = conditionMessage)
  expect_identical(bad("p99"), "the Hi-C colour max must be \"max\", \"p95\", \"p90\" or a number (it is p99).")
  # An old-style range: one message (a two-part one would stop the whole figure).
  expect_identical(bad(c(0, 50)), "the Hi-C colour max must be \"max\", \"p95\", \"p90\" or a number (it is 0, 50).")
})

# ---- On real data -------------------------------------------------------------------------------

test_that("on a real file, Rectangle reads and ranks the strip; p95, p90 and max differ (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- list(chrom = "chr2", start = 100000000, end = 105000000, gene = NA_character_)
  run <- function(cap, style = "rectangle") {
    s <- test_settings("mouse", list(zrange_cap = cap))
    e <- single_track("hic", test_rows("mouse")[["Ba"]], settings = s)
    e$hic_style <- style
    settings_hic(list(), e, s, 7.2, region)
  }
  out <- run("p95")
  m <- out$entry$data$matrix
  res <- out$entry$data$resolution
  # Read out to half a view width beyond each end (the track is half as high as wide).
  expect_lt(min(m$x), region$start - 2e6)
  expect_gt(max(m$y), region$end + 2e6)
  drawn <- .pt_hic_drawn_counts(m, "rectangle", region$start, region$end, res, out$entry$height / 7.2)
  drawn <- drawn[is.finite(drawn)]
  top <- function(o) o$dict[[.pt_shared_range_key("hic_zrange_", o$entry)]]
  expect_equal(top(out), unname(stats::quantile(drawn, 0.95)))
  expect_equal(top(run("p90")), unname(stats::quantile(drawn, 0.90)))
  expect_equal(top(run("max")), max(drawn))
  expect_lt(top(run("p90")), top(out))
  expect_lt(top(out), top(run("max")))
  # The view's square alone (what Auto used before): a higher 95th percentile, as the strip adds
  # long-range contacts, which are weaker.
  square <- top(run("p95", "triangle"))
  expect_lt(top(out), square)
})

test_that("the figure's colour key shows each choice's top (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  row <- test_rows("mouse")[["Ba"]]
  for (cap in list(NULL, "max", "p90")) {
    settings <- test_settings("mouse", list(zrange_cap = cap))
    f <- tempfile(fileext = ".svg")
    dict <- suppressMessages(suppressWarnings(render_all(list(single_track("hic", row, settings = settings)), region,
      settings, device_fn = function(w, h) svglite::svglite(f, width = w, height = h))))
    grDevices::dev.off()
    tops <- .pt_hic_colour_tops(dict)
    expect_length(tops, 1)
    expect_match(paste(readLines(f, warn = FALSE), collapse = "\n"),
                 paste0(">", format(signif(tops[[1]], 3), big.mark = ",", scientific = FALSE), "<"), fixed = TRUE,
                 info = format(cap))
    unlink(f)
  }
})

# ---- Settings round trips -----------------------------------------------------------------------

test_that("the Colour max box offers Auto, Max, 95th and 90th percentile, and keeps a saved number", {
  html <- function(cfg) as.character(.pt_display_settings_ui(cfg, "s_")$display)
  auto <- html(list())
  for (o in c("auto\" selected>Auto", "max\">Max", "p95\">95th percentile", "p90\">90th percentile"))
    expect_match(auto, paste0("<option value=\"", o, "</option>"), fixed = TRUE)
  expect_match(html(list(zrange_cap = "p90")), "<option value=\"p90\" selected>90th percentile", fixed = TRUE)
  expect_match(html(list(zrange_cap = 50L)), "<option value=\"50\" selected>50</option>", fixed = TRUE)
  expect_match(html(list(zrange_cap = 0.25)), "<option value=\"0.25\" selected>0.25</option>", fixed = TRUE)
  # A number can be typed in, nothing else.
  expect_match(auto, "\"createFilter\":\"^[0-9]*[.]?[0-9]+$\"", fixed = TRUE)
  expect_match(auto, "Auto is the 95th percentile on Linear and Max on Log", fixed = TRUE)
})

test_that("the settings panel's Colour max is read back as the setting", {
  read <- function(v, cfg = list()) .pt_collect_display_settings(list(s_zrange_cap = v), "s_", cfg)$zrange_cap
  expect_null(read("auto", list(zrange_cap = 50)))
  expect_identical(read("max"), "max")
  expect_identical(read("p95"), "p95")
  expect_identical(read("p90"), "p90")
  expect_identical(read("12.5"), 12.5)
  expect_null(read(""))
  expect_null(read(NA))
  # Not rendered yet: the saved one.
  expect_identical(.pt_collect_display_settings(list(), "s_", list(zrange_cap = "p90"))$zrange_cap, "p90")
})

test_that("project_settings.json keeps each Colour max; an old number or null still reads", {
  json <- tempfile(fileext = ".json")
  on.exit(unlink(json))
  read <- function(txt) { writeLines(txt, json); .pt_read_settings(json, .pt_default_settings()) }
  for (v in c("max", "p95", "p90")) {
    cfg <- expect_no_warning(read(sprintf('{"zrange_cap": "%s"}', v)))
    expect_identical(cfg$zrange_cap, v)
    expect_identical(.pt_settings_base(.pt_render_cfg(cfg, list(), tempdir()))$zrange_cap, v)
  }
  expect_identical(read('{"zrange_cap": 50}')$zrange_cap, 50L)
  expect_identical(.pt_settings_base(.pt_render_cfg(read('{"zrange_cap": 50}'), list(), tempdir()))$zrange_cap, 50L)
  expect_null(read('{"zrange_cap": null}')$zrange_cap)
  expect_null(.pt_settings_base(.pt_render_cfg(read('{"zrange_cap": null}'), list(), tempdir()))$zrange_cap)
  expect_warning(cfg <- read('{"zrange_cap": "p99"}'), "for zrange_cap; using the default")
  expect_null(cfg$zrange_cap)
})

test_that("Save writes the choice to project_settings.json and project_samples.txt, and both read back", {
  proj <- tempfile("pt-w95-")
  dir.create(file.path(proj, "data"), recursive = TRUE)
  on.exit(unlink(proj, recursive = TRUE))
  json <- file.path(proj, "project_settings.json")
  for (v in list("p90", "max", 50)) {
    # The wizard's Save: the collected display settings (no panel inputs) over the project's.
    settings <- c(list(data_dir = file.path(proj, "data"), samples = list(), gtf_path = ""),
                  .pt_collect_display_settings(list(), "", list(zrange_cap = v)))
    expect_true(.pt_write_project_files(settings, json))
    expect_equal(.pt_read_settings(json, .pt_default_settings())$zrange_cap, v)
    txt <- readLines(file.path(proj, "project_samples.txt"))
    expect_true(paste0("# zrange_cap: ", v) %in% txt, info = format(v))
    res <- import_samples_txt(file.path(proj, "project_samples.txt"), defaults = .pt_default_settings())
    expect_equal(.pt_coerce_auto_overrides(res$overrides, .pt_default_settings())$overrides$zrange_cap, v)
  }
})

test_that("a samples.txt Colour max: the choices and numbers import, anything else is named and ignored", {
  co <- function(v) .pt_coerce_auto_overrides(list(zrange_cap = v), .pt_default_settings())
  expect_identical(co("p95")$overrides$zrange_cap, "p95")
  expect_identical(co("50")$overrides$zrange_cap, 50)
  expect_null(co("auto")$overrides$zrange_cap)
  bad <- co("p99")
  expect_null(bad$overrides$zrange_cap)
  expect_identical(bad$ignored, "Ignoring '# zrange_cap: p99' (not max, p95, p90 or a number); using auto.")
  # Another auto setting still says "not a number".
  expect_match(.pt_coerce_auto_overrides(list(hic_height = "tall"), .pt_default_settings())$ignored,
               "(not a number)", fixed = TRUE)
})

# ---- The downloaded script ------------------------------------------------------------------------

hic_script <- function(cap, tops = NULL, scale = "linear", entries = NULL) {
  entries <- entries %||% list(list(viz = "hic", sample = "S", colour = "#000000", name = "HiC", path = "a.hic",
                                    hic_style = "rectangle"))
  assemble_download_script(entries, list(chrom = "chr2", start = 1000L, end = 5000L, gene = NA_character_),
                           list(orgdb_pkg = "org.Mm.eg.db", genome_label = "mm10", hic_norm = "KR",
                                hic_colour_scale = scale, zrange_cap = cap),
                           hic_tops = tops)
}

test_that("the script shows the Colour max chosen and the top it gave", {
  s <- hic_script("p95", c(hic_zrange_HiC = 41.5))
  expect_match(s, "  # In this figure: the 95th percentile, 41.5.\n  zrange_cap = \"p95\"", fixed = TRUE)
  expect_match(s, "# \"max\", \"p95\" or \"p90\" of the contacts drawn, or a number.", fixed = TRUE)
  back <- new.env()
  for (e in parse(text = s, keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("settings"))) eval(e, back)
  expect_identical(back$settings$zrange_cap, "p95")
  expect_match(hic_script(NULL, c(hic_zrange_HiC = 41.5)), "# In this figure: Auto, the 95th percentile, 41.5.\n  zrange_cap = NULL",
               fixed = TRUE)
  expect_match(hic_script(NULL, c(hic_zrange_HiC = 1234), "log"), "# In this figure: Auto, the maximum, 1,234.", fixed = TRUE)
  expect_match(hic_script("p90", c(hic_zrange_HiC = 0.0001235)), "the 90th percentile, 0.0001235.", fixed = TRUE)
  # A number is its own top; tops not known (none drawn): the choices only.
  for (s in list(hic_script(50, c(hic_zrange_HiC = 41.5)), hic_script("max"))) {
    expect_false(grepl("In this figure", s, fixed = TRUE))
    expect_match(s, "# Hi-C Colour max: NULL (Auto", fixed = TRUE)
  }
  expect_match(hic_script(50), "\n  zrange_cap = 50,\n", fixed = TRUE)
  # Two Hi-C tracks of different names: each named.
  two <- list(list(viz = "hic", sample = "S", colour = "#000000", name = "HiC", path = "a.hic"),
              list(viz = "hic", sample = "S", colour = "#000000", name = "loops", path = "b.hic"))
  expect_match(hic_script("max", c(hic_zrange_HiC = 41.5, hic_zrange_loops = 7), entries = two),
               "the maximum, 41.5 (HiC), 7 (loops).", fixed = TRUE)
  # No Hi-C track: no note.
  bw <- assemble_download_script(list(list(viz = "bw", sample = "S", colour = "#000000", name = "atac", path = "a.bw")),
                                 list(chrom = "chr2", start = 1000L, end = 5000L, gene = NA_character_),
                                 list(orgdb_pkg = "org.Mm.eg.db", genome_label = "mm10", zrange_cap = NULL))
  expect_false(grepl("# Hi-C Colour max:", bw, fixed = TRUE))
})

test_that("pt_script() names the top its Hi-C track draws (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  settings <- test_settings("mouse", list(zrange_cap = "p90"))
  fig <- structure(list(entries = list(single_track("hic", test_rows("mouse")[["Ba"]], settings = settings)),
                        region_params = test_region("mouse"), settings = settings), class = "pt_figure")
  tops <- suppressMessages(suppressWarnings(.pt_hic_figure_tops(fig)))
  expect_length(tops, 1)
  expect_gt(tops[[1]], 0)
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f))
  suppressMessages(suppressWarnings(pt_script(fig, f)))
  expect_match(paste(readLines(f), collapse = "\n"),
               paste0("# In this figure: the 90th percentile, ", format(tops[[1]], scientific = FALSE, big.mark = ","), "."),
               fixed = TRUE)
})

# ---- Normalisation errors say where to fix them --------------------------------------------------

test_that("a missing normalisation names the setting, its place and the values to pick, in 140 characters", {
  m <- .pt_hic_norm_missing("KR", 50000, c("VC", "SCALE", "NONE", "KR"))
  expect_identical(m, "this file has no KR normalisation at 50,000 bp bins: in Figure settings, Hi-C, set Normalisation (hic_norm) to VC, SCALE or NONE.")
  # No list of the file's own: NONE, which every file has.
  expect_identical(.pt_hic_norm_missing("KR", NULL, character(0)),
                   "this file has no KR normalisation: in Figure settings, Hi-C, set Normalisation (hic_norm) to NONE.")
  # A long list is cut to what fits, NONE always last.
  long <- .pt_hic_norm_missing("BALANCE", 1000000, c("SCALE", "VC_SQRT", "VC", "KR", "NONE", "INTER_KR", "GW_KR"))
  expect_lte(nchar(long), 140)
  expect_match(long, " or NONE.$")
})

test_that("a view with no values names the settings, their place and the values to set, in 140 characters", {
  expect_identical(.pt_hic_no_values("KR", 1000, 5000),
                   "no KR values at 1,000 bp bins in view: in Figure settings, Hi-C, set Normalisation (hic_norm) to NONE, or Resolution Fixed, Bin size 5000.")
  expect_identical(.pt_hic_no_values("KR", 1000, NULL),
                   "no KR values at 1,000 bp bins in view: in Figure settings, Hi-C, set Normalisation (hic_norm) to NONE.")
  expect_identical(.pt_hic_no_values("NONE", 1000, 5000),
                   "no contacts at 1,000 bp bins in view: in Figure settings, Hi-C, set Resolution Fixed, Bin size 5000.")
  expect_identical(.pt_hic_no_values("NONE", 1000, NULL), "no contacts at 1,000 bp bins in view.")
  # Too long with both: NONE only.
  expect_identical(.pt_hic_no_values("VC_SQRT", 1000000, 2500000),
                   "no VC_SQRT values at 1,000,000 bp bins in view: in Figure settings, Hi-C, set Normalisation (hic_norm) to NONE.")
})

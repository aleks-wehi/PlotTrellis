# Adding a track type the way a contributor does: a module (info_/settings_/render_
# functions) and nothing else, no list edited by hand. Two toy modules, defined here: a
# per-sample one (a BED file's features counted per bin, drawn as a signal; word for word the
# example module of the "Adding a track type" how-to) and a comparison (the log2 ratio of two
# bigWig tracks). Discovery is pointed at them by mocking the two functions every lookup
# goes through (.pt_viz_keys() and .pt_module_fn()), so everything that reads the
# declarations sees them, render_all() included, as for a module in R/.

info_bed_density <- function() list(
  kind = "sample", data_type = "bed", order = 35, label = "{name} density", filter = c(bed = "BED"),
  scale_bar = TRUE
)

settings_bed_density <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- settings$atac_height
  .pt_check_file(entry$path)  # a missing file fails the track
  lines <- tryCatch(readLines(entry$path), error = function(e) .pt_read_error(entry$path, e))
  lines <- lines[!grepl("^(#|track|browser)", lines)]  # header and comment lines
  raw <- if (length(lines) > 0) as.data.frame(data.table::fread(text = lines, header = FALSE, select = 1:3))
  chrom <- if (is.null(raw)) region_params$chrom else resolve_chrom(region_params$chrom, unique(as.character(raw[[1]])))
  if (is.null(chrom)) {
    entry$data <- list(no_chrom = TRUE)
    return(list(dict = dict, entry = entry))
  }
  if (is.null(raw)) raw <- data.frame(V1 = character(0), V2 = numeric(0), V3 = numeric(0))  # no features at all
  hits <- raw[raw[[1]] == chrom & raw[[3]] > region_params$start & raw[[2]] < region_params$end, ]
  breaks <- round(seq(region_params$start, region_params$end, length.out = 101))
  counts <- tabulate(findInterval((hits[[2]] + hits[[3]]) / 2, breaks, all.inside = TRUE), nbins = 100)
  entry$data <- list(no_chrom = FALSE, chrom = chrom,
                     gr = GRanges(chrom, IRanges(head(breaks, -1) + 1, tail(breaks, -1)), score = counts))
  # Every bed_density track of the same data name shares one value range.
  key <- .pt_shared_range_key("bed_density_range_", entry)
  dict[[key]] <- max(dict[[key]] %||% 0, counts)
  list(dict = dict, entry = entry)
}

render_bed_density <- function(entry, dict, ctx) {
  if (isTRUE(entry$data$no_chrom)) {
    .pt_status_text(ctx, "Chromosome not found")
    if (ctx$is_first) draw_scale_bar(ctx)
    return(invisible(NULL))
  }
  top <- dict[[.pt_shared_range_key("bed_density_range_", entry)]]
  range <- if (isTRUE(top > 0)) c(0, top * 10 / 9) else c(0, 1)
  plotSignal(data = entry$data$gr, params = ctx$region,
             chrom = entry$data$chrom, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
             x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height,
             fill = entry$colour, linecolor = entry$colour, range = range)
  .pt_heading(ctx, entry, paste0(entry$name, " density"))
  if (isTRUE(ctx$show_track_axes)) {
    draw_value_axis(range, ctx$track_x - 0.08, ctx$y, ctx$height, ctx$fontsize - 2, ctx$font, ctx$fontcol,
                    integer_ticks = TRUE)
  }
  if (ctx$is_first) draw_scale_bar(ctx)
}

info_bw_ratio <- function() list(
  kind = "comparison", data_type = "bw", order = 15, label = "{name} log2 ratio", api = "{name}_ratio",
  filter = c(comparison = "Comparisons")
)

settings_bw_ratio <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- settings$atac_height
  list(dict = dict, entry = entry)
}

render_bw_ratio <- function(entry, dict, ctx) {
  ctrl  <- load_bw_file(entry$control_path,   ctx$region_params)
  treat <- load_bw_file(entry$treatment_path, ctx$region_params)
  if (ctrl$no_chrom || treat$no_chrom || !identical(ctrl$chrom, treat$chrom)) {
    .pt_status_text(ctx, "Chromosome not found")
    return(invisible(NULL))
  }
  rp <- ctx$region_params
  breaks <- round(seq(rp$start, rp$end, length.out = 201))
  bins <- GRanges(ctrl$chrom, IRanges(head(breaks, -1) + 1, tail(breaks, -1)))
  mean_in <- function(gr) {
    cov <- GenomicRanges::coverage(gr, weight = "score")[[ctrl$chrom]]
    vapply(seq_along(bins), function(i) mean(as.numeric(cov[GenomicRanges::start(bins)[i]:GenomicRanges::end(bins)[i]])), 0)
  }
  lfc <- log2((mean_in(treat$gr) + 1) / (mean_in(ctrl$gr) + 1))
  lim <- max(abs(lfc), 0.1)
  range <- c(-lim, lim)
  up <- bins; up$score <- pmax(lfc, 0)
  dn <- bins; dn$score <- pmin(lfc, 0)
  plotSignal(data = up, params = ctx$region, chrom = ctrl$chrom, chromstart = rp$start, chromend = rp$end,
             x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height,
             fill = ctx$diff_colour_hi, linecolor = ctx$diff_colour_hi, range = range)
  plotSignal(data = dn, params = ctx$region, chrom = ctrl$chrom, chromstart = rp$start, chromend = rp$end,
             x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height,
             fill = ctx$diff_colour_lo, linecolor = ctx$diff_colour_lo, range = range)
  .pt_heading_pair(ctx, entry, paste0(entry$name, " log2 ratio"))
  if (isTRUE(ctx$show_track_axes)) {
    draw_value_axis(range, ctx$track_x - 0.08, ctx$y, ctx$height, ctx$fontsize - 2, ctx$font, ctx$fontcol)
  }
  if (ctx$is_first) draw_scale_bar(ctx)
}

TOYS <- c("bed_density", "bw_ratio")

# Makes the toy modules discoverable for the calling test. `seen` collects each toy
# render's ctx$is_first.
with_toys <- function(seen = new.env(), env = parent.frame()) {
  fns <- mget(paste0(rep(c("info_", "settings_", "render_"), each = 2), TOYS), envir = environment(with_toys))
  for (k in TOYS) local({
    key <- k
    draw <- fns[[paste0("render_", key)]]
    fns[[paste0("render_", key)]] <<- function(entry, dict, ctx) {
      seen[[key]] <- c(seen[[key]], ctx$is_first)
      draw(entry, dict, ctx)
    }
  })
  real_fn <- .pt_module_fn
  real_keys <- .pt_viz_keys
  local_mocked_bindings(
    .pt_viz_keys = function() c(real_keys(), TOYS),
    .pt_module_fn = function(prefix, key) {
      fn <- fns[[paste0(prefix, "_", key)]]
      if (is.null(fn)) real_fn(prefix, key) else fn
    },
    .env = env)
  seen
}

test_that("a new track type is offered, named, built into an entry and filed under its filter", {
  with_toys()
  rows <- list(
    list(display = "A", data = list(list(type = "bw", name = "atac", path = "A_atac.bw"),
                                    list(type = "bed", name = "peaks", path = "A_peaks.bed"))),
    list(display = "B", data = list(list(type = "bw", name = "atac", path = "B_atac.bw"))))
  items <- c(build_track_items(rows), build_comparison_track_items(rows))

  # Offered for a sample with its data type, not for one without.
  expect_true("bed_density_peaks_A" %in% names(items))
  expect_false(any(grepl("^bed_density_.*_B$", names(items))))
  expect_true(all(c("bw_ratio_atac_A_vs_B", "bw_ratio_atac_B_vs_A") %in% names(items)))
  # Detection order: after bed, and the ratio after the stacked diff.
  of_types <- function(types) names(Filter(function(it) it$type %in% types, items))
  expect_identical(of_types(c("bw", "bed", "bed_density")),
                   c("bw_atac_A", "bed_peaks_A", "bed_density_peaks_A", "bw_atac_B"))
  expect_identical(of_types(c("bw_diff_stacked", "bw_ratio")),
                   c("bw_diff_stacked_atac_A_vs_B", "bw_ratio_atac_A_vs_B", "bw_diff_stacked_atac_B_vs_A",
                     "bw_ratio_atac_B_vs_A"))

  # A label and an API name of its own.
  expect_identical(items$bed_density_peaks_A$label, "peaks density \u00B7 A")
  expect_identical(items$bw_ratio_atac_A_vs_B$label, "atac log2 ratio \u00B7 B (Cmp) vs A (Ref)")
  expect_false(anyDuplicated(vapply(items, `[[`, "", "label")) > 0)
  tab <- .pt_track_table(items)
  expect_identical(tab$track[tab$key == "bed_density_peaks_A"], "peaks_bed_density")
  expect_identical(tab$track[tab$key == "bw_ratio_atac_A_vs_B"], "atac_ratio")
  expect_false(anyDuplicated(paste(tab$track, tab$sample)) > 0)
  expect_identical(.pt_resolve_tracks(tab, c("peaks_bed_density@A", "atac_ratio@B vs A"), NULL),
                   c("bed_density_peaks_A", "bw_ratio_atac_A_vs_B"))

  # Entries with the sample, colour and paths.
  entries <- .pt_build_entries(c("bed_density_peaks_A", "bw_ratio_atac_A_vs_B"), items, c(A = "#E41A1C"),
                               list(reference_colour = "#111111", comparison_colour = "#222222"))
  expect_identical(entries[[1]][c("viz", "sample", "colour", "name", "path")],
                   list(viz = "bed_density", sample = "A", colour = "#E41A1C", name = "peaks", path = "A_peaks.bed"))
  expect_identical(entries[[2]][c("viz", "control_path", "treatment_path", "control_colour", "treatment_colour")],
                   list(viz = "bw_ratio", control_path = "A_atac.bw", treatment_path = "B_atac.bw",
                        control_colour = "#111111", treatment_colour = "#222222"))

  # Under its filter button, which already exists; no new button.
  filters <- .pt_type_filters()
  expect_identical(filters$bed_density, "bed")
  expect_identical(filters$bw_ratio, "comparison")
  infos <- .pt_viz_infos()
  expect_identical(as.character(tagList(.pt_filter_buttons())),
                   as.character(tagList(.pt_filter_buttons(infos[setdiff(names(infos), TOYS)]))))
  # The fixture helpers build it too.
  expect_identical(build_entry("bed_density", rows[[1]])$path, "A_peaks.bed")
  expect_identical(build_comparison_entry("bw_ratio", rows[[1]], rows[[2]])$treatment_path, "B_atac.bw")
})

test_that("a new track type renders through render_all(), with the scale bar as declared", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  seen <- with_toys()
  rows <- test_rows("human")
  items <- c(build_track_items(rows), build_comparison_track_items(rows))
  entries <- .pt_build_entries(c("bed_density_peaks_GM12878", "bw_ratio_atac_GM12878_vs_K562"), items,
                               .pt_sample_colours(rows), list())
  draw <- function(entries) {
    svg <- tempfile(fileext = ".svg")
    msgs <- character(0)
    # The ratio's lower layer is negative, which plotSignal() warns about.
    withCallingHandlers(
      suppressWarnings(render_all(entries, test_region("human"), test_settings("human"),
                                  device_fn = function(w, h) svglite::svglite(svg, width = w, height = h))),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    grDevices::dev.off()
    expect_false(any(grepl("track failed", msgs)))
    paste(readLines(svg, warn = FALSE), collapse = "\n")
  }
  svg <- draw(entries)
  expect_false(grepl("track failed", svg, fixed = TRUE))
  expect_match(svg, "peaks density", fixed = TRUE)
  expect_match(svg, "atac log2 ratio", fixed = TRUE)
  # bed_density declares scale_bar and is first; bw_ratio does not take it even when first.
  expect_identical(seen$bed_density, TRUE)
  expect_identical(seen$bw_ratio, FALSE)
  draw(rev(entries))
  expect_identical(seen$bed_density, c(TRUE, FALSE))
  expect_identical(seen$bw_ratio, c(FALSE, FALSE))
})

test_that("the example module counts every feature of a BED file with header and comment lines, and fails on a missing file", {
  bed <- tempfile(fileext = ".bed")
  on.exit(unlink(bed), add = TRUE)
  features <- c("chr1\t110\t120\ta", "chr1\t130\t140\tb", "chr1\t150\t160\tc")
  load <- function(lines = NULL, path = bed) {
    if (!is.null(lines)) writeLines(lines, bed)
    settings_bed_density(list(), list(name = "peaks", path = path), list(atac_height = 0.5), 6,
                         list(chrom = "chr1", start = 100L, end = 200L))$entry$data
  }
  total <- function(data) sum(data$gr$score)
  expect_identical(total(load(features)), 3L)
  expect_identical(total(load(c("browser position chr1:1-200", "track name=peaks", "# made by hand", features))), 3L)
  expect_identical(total(load(c(features[1], "# a comment between features", features[2:3]))), 3L)
  # A file with no features is an empty track, drawn as usual; a missing file is an error.
  for (data in list(load(character(0)), load("track name=peaks"))) {
    expect_false(data$no_chrom)
    expect_identical(total(data), 0L)
  }
  expect_error(load(path = tempfile()), "file not found")
  expect_error(load(path = NA_character_), "no data file given")
  expect_true(load(sub("chr1", "chr9", features))$no_chrom)
})

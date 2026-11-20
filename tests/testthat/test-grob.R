# pt_grob(): a figure as a grid grob, at its own size, that composes with other plots and
# keeps the same live text (in the same fonts) as pt_save() draws.
# Each label with its font family; pt_save()'s tidy step adds a generic fallback
# (", sans-serif") that a plain svglite file has not, so it is left out.
svg_texts <- function(f) {
  s <- grep("<text", readLines(f, warn = FALSE), value = TRUE)
  fam <- sub(", (sans-serif|serif|monospace)$", "", sub(".*font-family: ([^;']*);.*", "\\1", s))
  paste(fam, sub(".*>([^<]*)</text>.*", "\\1", s))
}
svg_sizes <- function(f) sort(as.numeric(sub(".*font-size: ([0-9.]+)px.*", "\\1",
                                               grep("<text", readLines(f, warn = FALSE), value = TRUE))))
# The x of every title `label` (a centred text), in points from the page's left.
title_x <- function(f, label) as.numeric(sub(".* x='([0-9.]+)'.*", "\\1",
                                             grep(paste0(">", label, "</text>"), readLines(f), value = TRUE)))
svg_inches <- function(f) {
  h <- readLines(f, n = 2)[2]
  as.numeric(strsplit(sub(".*width='([0-9.]+)pt' height='([0-9.]+)pt'.*", "\\1 \\2", h), " ")[[1]]) / 72
}
# Every label of `single`, as often as there, is in `composite` (which may hold more).
expect_labels_kept <- function(composite, single) {
  have <- table(svg_texts(composite)); want <- table(svg_texts(single))
  missing <- names(want)[is.na(have[names(want)]) | have[names(want)] < want]
  expect_identical(missing, character(0))
}
grob_project <- function() {
  dir <- tempfile("pt-grob-")
  dir.create(dir)
  list(dir = dir, proj = suppressWarnings(pt_project(file.path(dir, "proj"), pt_example())))
}

test_that("pt_grob() is the figure at its own size; grid.draw() draws it and a pt_figure", {
  skip_if_not_installed("org.Mm.eg.db")
  p <- grob_project()
  fig <- pt_figure(p$proj, "Eya2", c("ATAC", "genes"), title = "Locus Eya2")
  single <- pt_save(fig, file.path(p$dir, "Eya2.svg"))
  g <- pt_grob(fig)
  expect_s3_class(g, "pt_grob")
  expect_s3_class(g, "grob")
  expect_equal(c(g$width, g$height), svg_inches(single), tolerance = 1e-3)
  expect_equal(grid::widthDetails(g), grid::unit(g$width, "in"))
  expect_equal(grid::heightDetails(g), grid::unit(g$height, "in"))

  # Drawn twice side by side: the labels of two figures, in their fonts, the second
  # figure one figure width to the right (not drawn into the first one's page).
  f <- file.path(p$dir, "two.svg")
  svglite::svglite(f, width = 2 * g$width, height = g$height + 1)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  for (j in 1:2) {
    grid::pushViewport(grid::viewport(layout.pos.col = j))
    grid::grid.draw(g)
    grid::upViewport()
  }
  grDevices::dev.off()
  expect_equal(sort(svg_texts(f)), sort(rep(svg_texts(single), 2)))
  expect_equal(diff(title_x(f, "Locus Eya2")[1:2]), 72 * g$width, tolerance = 1e-3)
  # Settings of a viewport around it do not change its text.
  f <- file.path(p$dir, "cex.svg")
  svglite::svglite(f, width = g$width, height = g$height)
  grid::pushViewport(grid::viewport(gp = grid::gpar(cex = 2, fontface = "italic")))
  grid::grid.draw(g)
  grDevices::dev.off()
  expect_equal(sort(svg_texts(f)), sort(svg_texts(single)))
  expect_equal(svg_sizes(f), svg_sizes(single))
  expect_false(any(grepl("font-style: italic", readLines(f))))
  # A parent alpha (which cannot be undone, alpha stays at most 1) and lex: drawn, same text.
  f <- file.path(p$dir, "alpha.svg")
  svglite::svglite(f, width = g$width, height = g$height)
  grid::pushViewport(grid::viewport(gp = grid::gpar(alpha = 0.5, lex = 2)))
  expect_no_error(grid::grid.draw(g))
  grDevices::dev.off()
  expect_equal(svg_sizes(f), svg_sizes(single))
  # R's pdf() device does not know the font: drawn in its own, one message, no viewport
  # left behind (see the every-track test below). On macOS without DejaVu Sans the
  # figure is in "sans", which pdf() draws.
  grDevices::pdf(NULL)
  pdf_dev <- grDevices::dev.cur()
  if (g$font %in% names(grDevices::pdfFonts())) {
    expect_no_message(grid::grid.draw(g))
  } else {
    expect_message(grid::grid.draw(g), "device = cairo_pdf")
  }
  expect_identical(grid::current.vpPath(), NULL)
  grDevices::dev.off(pdf_dev)
  # The pt_figure itself.
  f <- file.path(p$dir, "fig.svg")
  svglite::svglite(f, width = g$width, height = g$height)
  grid::grid.draw(fig)
  grDevices::dev.off()
  expect_equal(sort(svg_texts(f)), sort(svg_texts(single)))

  expect_error(pt_grob("Eya2"), "takes a pt_figure")
  local_mocked_bindings(.pt_capture = function(fig, ext, dpi) stop("boom"))
  expect_error(pt_grob(fig), "Could not draw the figure of Eya2: boom")
})

test_that("a pt_grob composes with ggplot2 in patchwork and cowplot, text kept", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("cowplot")
  p <- grob_project()
  fig <- pt_figure(p$proj, "Eya2", c("ATAC", "genes"), title = "Locus Eya2")
  single <- pt_save(fig, file.path(p$dir, "Eya2.svg"))
  g <- pt_grob(fig)
  scatter <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::ggtitle("a ggplot beside it")

  # The page about half an inch bigger than the figure, as ?pt_grob says.
  pw <- patchwork::wrap_elements(full = g) + scatter +
    patchwork::plot_layout(widths = grid::unit(c(g$width, 1), c("in", "null")))
  f <- file.path(p$dir, "patchwork.svg")
  expect_no_warning(ggplot2::ggsave(f, pw, width = g$width + 4.5, height = g$height + 0.5))
  expect_labels_kept(f, single)
  expect_true(any(grepl(">a ggplot beside it<", readLines(f))))
  expect_false(any(grepl("<image", readLines(f))))
  # Two copies side by side: the second one figure width or more to the right.
  f <- file.path(p$dir, "patchwork2.svg")
  ggplot2::ggsave(f, patchwork::wrap_elements(full = g) + patchwork::wrap_elements(full = g),
                  width = 2 * g$width + 0.5, height = g$height + 0.5)
  expect_gte(diff(title_x(f, "Locus Eya2")[1:2]), 72 * g$width)

  cw <- cowplot::plot_grid(g, scatter, labels = c("A", "B"), rel_widths = c(g$width, 4))
  f <- file.path(p$dir, "cowplot.svg")
  expect_no_warning(ggplot2::ggsave(f, cw, width = g$width + 4.5, height = g$height + 0.5))
  expect_labels_kept(f, single)
  expect_true(any(grepl(">B<", readLines(f))))
  # A plain ggsave() of a .pdf uses R's pdf(): the figure in its font, with no warning,
  # the message saying how to keep DejaVu Sans (unless pdf() knows the font, as "sans"
  # on macOS without DejaVu Sans), and the figure's text in Helvetica.
  f <- file.path(p$dir, "patchwork-default.pdf")
  save_default <- function() expect_no_warning(ggplot2::ggsave(f, pw, width = g$width + 4.5,
                                                              height = g$height + 0.5))
  if (g$font %in% names(grDevices::pdfFonts())) expect_no_message(save_default())
  else expect_message(save_default(), "device = cairo_pdf")
  pdf_text <- paste(rawToChar(readBin(f, "raw", file.size(f)), multiple = TRUE), collapse = "")
  expect_true(grepl("/BaseFont /Helvetica", pdf_text, fixed = TRUE, useBytes = TRUE))
  expect_false(grepl("DejaVu", pdf_text, fixed = TRUE, useBytes = TRUE))
})

# R's pdf() and postscript() devices, which lack the default DejaVu Sans, draw a pt_grob
# of every track type in their own font, with one message and no error or warning
# (postscript() warns only that it cannot draw semi-transparency, as for any plot).
test_that("a pt_grob of every track type draws on R's pdf() and postscript() devices", {
  skip_if_not_installed("org.Mm.eg.db")
  p <- grob_project()
  fig <- pt_figure(p$proj, "Eya2", unique(pt_tracks(p$proj)$track), comparison = c("LP", "Basal"))
  expect_setequal(vapply(fig$entries, `[[`, "", "viz"), c(
    .pt_viz_keys()[.pt_viz_keys() != "genome_axis"]))
  # A generic family ("sans", the default on macOS without DejaVu Sans) is in
  # postscriptFonts() but postscript() draws only its own family: drawn in that, no message.
  fig_sans <- fig
  fig_sans$settings$font <- "sans"
  g_sans <- pt_grob(fig_sans)
  f <- file.path(p$dir, "sans.ps")
  grDevices::postscript(f, width = g_sans$width, height = g_sans$height)
  expect_no_message(suppressWarnings(expect_no_error(grid::grid.draw(g_sans)), classes = "warning"))
  grDevices::dev.off()
  expect_identical(grid::current.vpPath(), NULL)
  g <- pt_grob(fig)
  expect_identical(g$font, .pt_device_font("DejaVu Sans"))
  skip_if(g$font %in% names(grDevices::pdfFonts()), "the figure's font is one pdf() knows")
  for (dev in c("pdf", "postscript")) {
    f <- file.path(p$dir, paste0("all.", if (dev == "pdf") "pdf" else "ps"))
    if (dev == "pdf") grDevices::pdf(f, width = g$width, height = g$height, compress = FALSE)
    else grDevices::postscript(f, width = g$width, height = g$height)
    msgs <- character(0)
    withCallingHandlers(
      expect_no_error(grid::grid.draw(g)),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
      },
      warning = function(w) {
        if (dev == "postscript" && grepl("semi-transparency", conditionMessage(w)))
          invokeRestart("muffleWarning")
      })
    grDevices::dev.off()
    expect_identical(grid::current.vpPath(), NULL)
    expect_length(msgs, 1)
    expect_match(msgs, paste0("R's ", dev, "\\(\\) device has no DejaVu Sans.*device = cairo_pdf"))
    txt <- paste(rawToChar(readBin(f, "raw", file.size(f)), multiple = TRUE), collapse = "")
    # Helvetica (a postscript file names it in its prolog anyway), the samples' labels
    # (short ones: pdf() splits longer text at kerning pairs), and no track drawn as failed.
    if (dev == "pdf") expect_true(grepl("/BaseFont /Helvetica", txt, fixed = TRUE, useBytes = TRUE))
    expect_false(grepl("DejaVu", txt, fixed = TRUE, useBytes = TRUE))
    for (label in c("(Basal)", "(LP)"))
      expect_true(grepl(label, txt, fixed = TRUE, useBytes = TRUE), info = paste(dev, label))
    expect_false(grepl("track failed", txt, fixed = TRUE, useBytes = TRUE))
  }
})

test_that("a pt_grob + a ggplot says to use patchwork, instead of a silent NULL (bug 266, no fixture data)", {
  skip_if_not_installed("ggplot2")
  g <- grid::grob(cl = "pt_grob")
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_error(g + p, "patchwork::wrap_elements(full = g) + p", fixed = TRUE)
  expect_error(g + 1, "patchwork::wrap_elements", fixed = TRUE)
  # p + g stays patchwork's, also through this method (ggplot2 before 4.0 dispatches here).
  skip_if_not_installed("patchwork")
  expect_s3_class(p + g, "patchwork")
  expect_s3_class(`+.pt_grob`(p, g), "patchwork")
})

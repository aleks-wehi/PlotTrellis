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
  fig <- pt_figure(p$proj, "Elf5", c("atac", "genes"), title = "Locus Elf5")
  single <- pt_save(fig, file.path(p$dir, "Elf5.svg"))
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
  expect_equal(diff(title_x(f, "Locus Elf5")[1:2]), 72 * g$width, tolerance = 1e-3)
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
  # R's pdf() device does not know the font: a clear stop, and no viewport left behind.
  grDevices::pdf(NULL)
  pdf_dev <- grDevices::dev.cur()
  expect_warning(expect_error(grid::grid.draw(g), "draw it on a cairo device"), "draw it on a cairo device")
  expect_identical(grid::current.vpPath(), NULL)
  grDevices::dev.off(pdf_dev)
  # The pt_figure itself.
  f <- file.path(p$dir, "fig.svg")
  svglite::svglite(f, width = g$width, height = g$height)
  grid::grid.draw(fig)
  grDevices::dev.off()
  expect_equal(sort(svg_texts(f)), sort(svg_texts(single)))

  expect_error(pt_grob("Elf5"), "takes a pt_figure")
  local_mocked_bindings(.pt_capture = function(fig, ext, dpi) stop("boom"))
  expect_error(pt_grob(fig), "Could not draw the figure of Elf5: boom")
})

test_that("a pt_grob composes with ggplot2 in patchwork and cowplot, text kept", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("cowplot")
  p <- grob_project()
  fig <- pt_figure(p$proj, "Elf5", c("atac", "genes"), title = "Locus Elf5")
  single <- pt_save(fig, file.path(p$dir, "Elf5.svg"))
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
  expect_gte(diff(title_x(f, "Locus Elf5")[1:2]), 72 * g$width)

  cw <- cowplot::plot_grid(g, scatter, labels = c("A", "B"), rel_widths = c(g$width, 4))
  f <- file.path(p$dir, "cowplot.svg")
  expect_no_warning(ggplot2::ggsave(f, cw, width = g$width + 4.5, height = g$height + 0.5))
  expect_labels_kept(f, single)
  expect_true(any(grepl(">B<", readLines(f))))
  # patchwork drops errors raised while it draws: the warning still says what to do.
  f <- file.path(p$dir, "patchwork-default.pdf")
  expect_warning(ggplot2::ggsave(f, pw, width = g$width + 4.5, height = g$height + 0.5, device = grDevices::pdf),
                 "draw it on a cairo device")
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

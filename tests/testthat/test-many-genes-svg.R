# The paper's many-genes chunk: pt_figure() + pt_save() per gene in a plain loop, then
# xml2 (an external package, Suggests) stitches the SVGs into one grid. Same code as the
# chunk handed to the owner (solutions/review/W89/chunk.R in the monorepo; change both
# together), only writing into a temporary folder, with a transcripts track so that the
# figures differ in height (Cxcl14 has fewer transcripts).
test_that("a loop of pt_save() SVGs stitches into one vector SVG with xml2", {
  skip_if_not_installed("xml2")
  skip_if_not_installed("org.Mm.eg.db")
  if (!"package:xml2" %in% search()) on.exit(detach("package:xml2"), add = TRUE)
  library(xml2)
  dir <- tempfile("pt-many-genes-")
  dir.create(dir)
  proj <- suppressWarnings(pt_project(file.path(dir, "my-project"), pt_example()))
  genes <- c("Eya2", "Acot7", "Cxcl14")
  svgs <- c()
  for (g in genes)
    svgs[g] <- pt_save(pt_figure(proj, g, c("ATAC", "transcripts_all", "genes"), title = g), file.path(dir, paste0(g, ".svg")))
  expect_true(all(file.exists(svgs)))   # every gene's own SVG is kept

  figs <- lapply(svgs, read_xml); ncol <- min(2, length(figs))
  wh <- sapply(figs, function(f) as.numeric(sub("pt", "", c(xml_attr(f, "width"), xml_attr(f, "height")))))
  w <- max(wh[1, ]); h <- max(wh[2, ])
  big <- read_xml(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%gpt" height="%gpt" viewBox="0 0 %1$g %2$g"/>',
                          ncol * w, ceiling(length(figs) / ncol) * h))
  for (i in seq_along(figs)) {
    xml_attr(figs[[i]], "x") <- (i - 1) %% ncol * w; xml_attr(figs[[i]], "y") <- (i - 1) %/% ncol * h
    xml_attr(figs[[i]], "width") <- wh[1, i];        xml_attr(figs[[i]], "height") <- wh[2, i]
    xml_add_child(big, figs[[i]])
  }
  out <- file.path(dir, "all_genes.svg")
  write_xml(big, out)

  # Valid XML, one page of 2 x 2 cells, sized in points like the single figures.
  x <- read_xml(out)
  expect_equal(xml_name(x), "svg")
  expect_equal(xml_attr(x, "width"), paste0(2 * w, "pt"))
  expect_equal(xml_attr(x, "height"), paste0(2 * h, "pt"))
  # Three nested figures, in grid order, each still vector (shapes and live text, no raster).
  kids <- xml_find_all(x, "/*[local-name()='svg']/*[local-name()='svg']")
  expect_length(kids, 3)
  expect_equal(as.numeric(xml_attr(kids, "x")), c(0, w, 0))
  expect_equal(as.numeric(xml_attr(kids, "y")), c(0, 0, h))
  # Each figure keeps its own size (top-left in its cell, not stretched), its viewBox and
  # its xlink namespace (for xlink:href, should a track ever embed an image).
  expect_gt(length(unique(wh[2, ])), 1)
  expect_equal(as.numeric(xml_attr(kids, "height")), unname(wh[2, ]))
  expect_equal(as.numeric(xml_attr(kids, "width")), unname(wh[1, ]))
  expect_true(all(!is.na(xml_attr(kids, "viewBox"))))
  expect_true("http://www.w3.org/1999/xlink" %in% xml_ns(x))
  for (i in seq_along(genes)) {
    txt <- xml_text(xml_find_all(kids[[i]], ".//*[local-name()='text']"))
    expect_true(genes[i] %in% txt)   # the figure's title
    expect_gt(length(xml_find_all(kids[[i]], ".//*[local-name()='polygon' or local-name()='rect' or local-name()='line']")), 50)
  }
  expect_length(xml_find_all(x, "//*[local-name()='image']"), 0)
  # Non-trivial: about the three files together, none dropped.
  expect_gt(file.size(out), 0.9 * sum(file.size(svgs)))
  # svglite's clip-path ids repeat across figures; they are named after the clip rectangle,
  # so a repeated id always means the same rectangle (in its own figure's coordinates).
  clips <- xml_find_all(x, "//*[local-name()='clipPath']")
  def <- vapply(clips, function(n) as.character(n), "")
  ids <- xml_attr(clips, "id")
  expect_true(all(tapply(def, ids, function(d) length(unique(d)) == 1)))
})

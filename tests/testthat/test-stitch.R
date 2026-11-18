# pt_stitch(): figures stacked top to bottom in one .svg, .pdf or .png, each as pt_save() draws it.
test_that("pt_stitch() stacks figures in one file, as pt_save() draws each", {
  skip_if_not_installed("org.Mm.eg.db")
  dir <- tempfile("pt-stitch-")
  dir.create(dir)
  proj <- suppressWarnings(pt_project(file.path(dir, "proj"), pt_example()))
  # Different heights (Cat has fewer transcripts); the genes track's labels take their font
  # from a viewport around plotGenes(), which a replay through grid.draw() lost.
  figs <- list()
  for (g in c("Elf5", "Cat"))
    figs[[g]] <- pt_figure(proj, g, c("atac", "transcripts_all", "genes"), title = g)
  singles <- vapply(names(figs), function(g) pt_save(figs[[g]], file.path(dir, paste0(g, ".svg"))), "")

  root_wh <- function(f) {
    head <- readLines(f, n = 2)[2]
    as.numeric(sub(".*width='([0-9.]+)pt' height='([0-9.]+)pt'.*", "\\1 \\2", head) |> strsplit(" ") |> unlist())
  }
  texts <- function(f) {
    s <- grep("<text", readLines(f, warn = FALSE), value = TRUE)
    sort(paste(sub(".*font-family: ([^;']*);.*", "\\1", s), sub(".*>([^<]*)</text>.*", "\\1", s)))
  }
  wh <- sapply(singles, root_wh)

  dev_before <- grDevices::dev.cur()
  out <- pt_stitch(figs, file.path(dir, "all.svg"))
  expect_equal(grDevices::dev.cur(), dev_before)          # no device left open or switched
  expect_equal(out, file.path(dir, "all.svg"))
  s <- readLines(out, warn = FALSE)
  expect_equal(sum(grepl("^<svg ", s)), 1)                  # one flat SVG, not nested files
  expect_equal(root_wh(out), c(max(wh[1, ]), sum(wh[2, ])), tolerance = 0.02)
  expect_false(any(grepl("<image", s)))                     # vector, no bitmap
  # The same labels in the same fonts as the two single figures (live text).
  expect_equal(texts(out), sort(c(texts(singles[1]), texts(singles[2]))))
  expect_true(any(grepl(">Elf5-201<", s)) && any(grepl(">Cat-202<", s)))
  # Cat below Elf5: its title is lower on the page than Elf5's height.
  y_of <- function(label) as.numeric(sub(".* y='([0-9.]+)'.*", "\\1", grep(paste0(">", label, "</text>"), s, value = TRUE)[1]))
  expect_lt(y_of("Elf5"), wh[2, 1])
  expect_gt(y_of("Cat"), wh[2, 1])

  # Side by side (ncol = 2): one row, as wide as both, as tall as the taller.
  side <- pt_stitch(figs, file.path(dir, "side.svg"), ncol = 2)
  expect_equal(root_wh(side), c(sum(wh[1, ]), max(wh[2, ])), tolerance = 0.02)
  expect_equal(texts(side), texts(out))
  expect_gt(as.numeric(sub(".* x='([0-9.]+)'.*", "\\1", grep(">Cat-202</text>", readLines(side), value = TRUE))), wh[1, 1])

  pdf <- pt_stitch(figs, file.path(dir, "all.pdf"))
  # cairo keeps the page objects in compressed streams: inflate every stream, then search.
  raw <- readBin(pdf, "raw", file.size(pdf))
  parts <- c(list(raw), Map(function(a, b) {
    a <- a + if (raw[a + 7] == as.raw(13)) 9 else 8
    tryCatch(memDecompress(raw[a:(b - 1)], "gzip"), error = function(e) raw(0))
  }, grepRaw("[^d]stream\r?\n", raw, all = TRUE), grepRaw("\r?\nendstream", raw, all = TRUE)))
  find <- function(re) unlist(lapply(parts, function(x) vapply(grepRaw(re, x, value = TRUE, all = TRUE), rawToChar, "")))
  expect_length(find("/Type */Page[^s]"), 1)                                         # one page
  box <- as.numeric(strsplit(trimws(gsub("[^0-9. ]", "", find("/MediaBox *\\[[^]]*\\]"))), " +")[[1]])
  expect_equal(box[3:4], c(max(wh[1, ]), sum(wh[2, ])), tolerance = 0.02)
  expect_gt(length(find("/FontFile[23]?")), 0)                                      # embedded fonts: live text

  png <- pt_stitch(figs, file.path(dir, "all.png"), dpi = 50)
  expect_equal(rawToChar(readBin(png, "raw", 4)[2:4]), "PNG")
  px <- readBin(png, "raw", 24)
  dims <- c(sum(as.integer(px[17:20]) * 256^(3:0)), sum(as.integer(px[21:24]) * 256^(3:0)))
  expect_equal(dims, round(c(max(wh[1, ]), sum(wh[2, ])) / 72 * 50), tolerance = 0.01)

  expect_error(pt_stitch(figs, file.path(dir, "all.jpg")), "writes .pdf, .svg or .png")

  # Plain errors, before any drawing where possible; a failing figure is named.
  expect_error(pt_stitch(list(), file.path(dir, "x.pdf")), "figs must be a list of pt_figure")
  expect_error(pt_stitch(list(figs[[1]], NULL), file.path(dir, "x.pdf")), "item 2 is not one")
  expect_error(pt_stitch(figs, file.path(dir, "no-such-folder", "x.pdf")), "does not exist")
  local_mocked_bindings(.pt_capture = function(fig, ext, dpi) {
    if (fig$gene == "Cat") stop("boom")
    structure(grid::gTree(), pt_size = c(9, 60))
  })
  expect_error(pt_stitch(figs, file.path(dir, "x.pdf")), "Figure 2 \\(Cat\\): boom")
  # Two 60-inch figures: past cairo's 32767-pixel PNG limit at 300 dpi, fine at 100.
  tall <- list(figs[[1]], figs[[1]])
  expect_error(pt_stitch(tall, file.path(dir, "tall.png")), "Too big for one PNG at 300 dpi")
  expect_no_error(pt_stitch(tall, file.path(dir, "tall.png"), dpi = 100))
  expect_equal(grDevices::dev.cur(), dev_before)
})

# Cell layout with made-up figure sizes (the capture is mocked, nothing is rendered): each
# figure's viewport is recorded as its top-left corner and size in inches on the page.
test_that("pt_stitch(ncol =) fills rows left to right, each figure top-left in its cell", {
  skip_if_not_installed("org.Mm.eg.db")
  dir <- tempfile("pt-stitch-ncol-")
  dir.create(dir)
  proj <- suppressWarnings(pt_project(file.path(dir, "proj"), pt_example()))
  fig <- pt_figure(proj, "Elf5", "genes")
  figs <- list(fig, fig, fig)
  sizes <- list(c(4, 5), c(3, 2), c(6, 1))   # widths and heights in inches, in order
  drawn <- list()
  local_mocked_bindings(.pt_capture = function(fig, ext, dpi) {
    size <- sizes[[length(drawn) %% length(sizes) + 1]]
    drawn[[length(drawn) + 1]] <<- size
    record <- function() {
      tl <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
      wh <- c(grid::convertWidth(grid::unit(1, "npc"), "in", TRUE),
              grid::convertHeight(grid::unit(1, "npc"), "in", TRUE))
      at[[length(at) + 1]] <<- round(c(tl$x, tl$y, wh), 3)
    }
    structure(list(children = list(list(expr = quote(record()), list = list(record = record)))),
              pt_size = size)
  })
  page_in <- function(f) {
    h <- readLines(f, n = 2)[2]
    as.numeric(strsplit(sub(".*width='([0-9.]+)pt' height='([0-9.]+)pt'.*", "\\1 \\2", h), " ")[[1]]) / 72
  }

  # ncol = 2, 3 figures: columns 6 (max of 4 and 6) and 3 wide, rows 5 (max of 5 and 2) and 1
  # tall, so the page is 9 x 6 in; the last cell (right of the third figure) stays empty.
  at <- list(); drawn <- list()
  f <- pt_stitch(figs, file.path(dir, "grid.svg"), ncol = 2)
  expect_equal(page_in(f), c(9, 6), tolerance = 1e-3)
  expect_equal(at, list(c(0, 6, 4, 5),     # row 1, column 1
                        c(6, 6, 3, 2),     # row 1, column 2, starts after column 1's width
                        c(0, 1, 6, 1)))    # row 2, column 1, below row 1's tallest figure

  # ncol larger than the number of figures: one row.
  at <- list(); drawn <- list()
  f <- pt_stitch(figs, file.path(dir, "row.svg"), ncol = 5)
  expect_equal(page_in(f), c(13, 5), tolerance = 1e-3)
  expect_equal(at, list(c(0, 5, 4, 5), c(4, 5, 3, 2), c(7, 5, 6, 1)))
  # The default is one column, as before.
  at <- list(); drawn <- list()
  f <- pt_stitch(figs, file.path(dir, "col.svg"))
  expect_equal(page_in(f), c(6, 8), tolerance = 1e-3)
  expect_equal(at, list(c(0, 8, 4, 5), c(0, 3, 3, 2), c(0, 1, 6, 1)))

  for (bad in list(0, -1, 1.5, NA, NA_real_, "2", c(1, 2), NULL))
    expect_error(pt_stitch(figs, file.path(dir, "bad.svg"), ncol = bad), "ncol must be a whole number")

  # The PNG limit uses the whole page, wide as well as tall: three 40 x 1 in figures fit
  # stacked (40 x 3 in, 12000 x 900 px at 300 dpi) but not side by side (120 in, 36000 px).
  sizes <- list(c(40, 1))
  expect_error(pt_stitch(figs, file.path(dir, "wide.png"), ncol = 3), "Too big for one PNG at 300 dpi")
  expect_no_error(pt_stitch(figs, file.path(dir, "wide.png"), ncol = 1))
})

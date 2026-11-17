# In-track status texts that modules draw themselves (no retained or annotated transcripts, the
# protein-domain statuses) stay inside a narrow track: drawn through no_data_text(), which fits them to the
# space below the track's heading (wrapped, smaller if need be).

# The left and right ends (inches) of the centred texts a track draws, its heading left out.
status_text_ends <- function(draw, track_x = 0.5, track_w = 1.5, height = 1, heading) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  page_w <- track_x + track_w + 0.5
  svglite::svglite(f, width = page_w, height = height + 0.5)
  pageCreate(width = page_w, height = height + 0.5, default.units = "inches", showGuides = FALSE,
             xgrid = 0, ygrid = 0)
  draw(list(track_x = track_x, track_w = track_w, y = 0.25, height = height, fontsize = 8, font = "sans",
            fontcol = "black", region_params = list(gene = "Gene1")))
  grDevices::dev.off()
  txt <- grep("^<text ", readLines(f, warn = FALSE), value = TRUE)
  # Status texts are centred in the track; the heading and a subtitle are left-aligned.
  txt <- txt[!grepl(heading, txt, fixed = TRUE) & grepl("text-anchor='middle'", txt, fixed = TRUE)]
  x <- as.numeric(sub("^<text x='([0-9.]+)'.*$", "\\1", txt)) / 72
  len <- as.numeric(sub("^.* textLength='([0-9.]+)px'.*$", "\\1", txt)) / 72
  list(n = length(txt), left = min(x - len / 2), right = max(x + len / 2))
}

expect_inside_track <- function(ends, track_x = 0.5, track_w = 1.5) {
  expect_gt(ends$n, 0)
  expect_gte(ends$left, track_x - 1e-3)
  expect_lte(ends$right, track_x + track_w + 1e-3)
}

test_that("'No retained/annotated transcripts' stays inside a narrow transcripts track", {
  expect_inside_track(status_text_ends(function(ctx) render_transcripts_all(list(data = list(tx = data.frame())), NULL, ctx),
                                       heading = ">Transcripts<"))
  entry <- list(data = list(n_excluded = 0, kept = data.frame()))
  expect_inside_track(status_text_ends(function(ctx) render_transcripts(entry, NULL, ctx),
                                       heading = ">Transcripts (filtered)<"))
  expect_inside_track(status_text_ends(function(ctx) render_transcripts_annotated(entry, NULL, ctx),
                                       heading = ">Transcripts annotated (filtered)<"))
})

test_that("protein-domain status texts stay inside a narrow track", {
  heading <- ">Canonical protein annotation<"
  local_mocked_bindings(gene_protein_domains = function(...) list(status = "no_canonical_transcript"))
  expect_inside_track(status_text_ends(function(ctx) render_protein_domains(list(), NULL, ctx), heading = heading))
  local_mocked_bindings(gene_protein_domains = function(...) list(
    status = "ok", canonical_tx = "Gene1-201", method = "ccds_id_join",
    domains = data.frame(domain_id = character(0), start = numeric(0))))
  ends <- status_text_ends(function(ctx) render_protein_domains(list(), NULL, ctx), heading = heading)
  expect_inside_track(ends)
})

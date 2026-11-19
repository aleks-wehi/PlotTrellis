# Effect tests for reference_colour / comparison_colour, whose only consumer is the app's
# server, not render_all(): they are resolved onto comparison entries by .pt_build_entries().
# test-composer-settings.R checks every setting reaches the render; this file checks these
# two change what it draws. uniprot_reviewed_only's effect is test-main-api.R's (the Gene
# info UniProt link) and test-main-server.R's (the live value reaches it); gene_flank_pct's
# is test-composer-settings.R's. Drives the real server (pattern of test-main-server.R).
# Needs the mouse fixture.

effects_composer <- function() {
  proj <- tempfile("pt-effects-"); dir.create(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  list(proj = proj, app = suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()), experienced = TRUE)))
}

test_that("reference_colour/comparison_colour colour the comparison tracks (mouse dataset)", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- effects_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- recording_session()
  # Ba (reference) against LP (comparison): the track bw_diff_stacked_atac_Ba_vs_LP
  # (R/viz-track-detection.R); the page sends the comparison first.
  args <- list(chr = "chr2", s = 103400000, e = 103460000, picked = list("LP", "Ba"),
               tracks = list(list(type = "atac_diff")))
  svg_of <- function(id) {
    session <- rec$session
    session$setInputs(pt_rpc = as.character(jsonlite::toJSON(list(list(id = id, op = "render", args = args)),
                                                             auto_unbox = TRUE, digits = NA)))
    r <- Filter(function(m) identical(m$type, "pt_rpc_reply") && identical(m$message$id, id), rec$custom)
    fig <- jsonlite::fromJSON(r[[length(r)]]$message$json, simplifyVector = FALSE)
    expect_null(fig$error)
    fig$svg
  }
  shiny::testServer(cp$app, {
    session$flushReact()
    def <- svg_of("r1")
    expect_true(grepl("#4D4D4D", def, fixed = TRUE))
    expect_true(grepl("#56B4E9", def, fixed = TRUE))

    session$setInputs(s_reference_colour = "#123456", s_comparison_colour = "#ABCDEF")
    alt <- svg_of("r2")
    expect_true(grepl("#123456", alt, fixed = TRUE))
    expect_true(grepl("#ABCDEF", alt, fixed = TRUE))
    expect_false(grepl("#4D4D4D", alt, fixed = TRUE))
    expect_false(grepl("#56B4E9", alt, fixed = TRUE))
  }, session = rec$session)
})

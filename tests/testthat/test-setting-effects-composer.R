# Phase I.3 step 28 (stage B): effect tests for the settings whose only
# consumer is the composer server, not render_all(): reference_colour /
# comparison_colour (resolved onto comparison entries by build_entries()) and
# uniprot_reviewed_only (the Gene info panel's UniProt link-out).
# test-composer-settings.R checks these values reach the composer; this file
# checks they change what it produces. gene_flank_pct's effect is covered
# there ("a composer flank change re-resolves the current gene's region").
# Needs the mouse fixture.

effects_composer <- function() {
  proj <- tempfile("pt-effects-"); dir.create(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  list(proj = proj, app = suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))
}

test_that("reference_colour/comparison_colour colour the comparison tracks (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- effects_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    # Comparison item key format: R/track_detection.R (<viz>_<name>_<ref>_vs_<cmp>).
    key <- "bw_diff_stacked_atac_Ba_vs_LP"
    session$setInputs(view_mode = "compose", region_mode = "coords", region_chr = "chr2",
                      region_range = c(103400000L, 103460000L), figure_tracks = key)
    svg_now <- function() paste(readLines(state$svg_path, warn = FALSE), collapse = "\n")
    session$setInputs(generate = 1)
    expect_null(state$plot_error)
    def <- svg_now()
    expect_true(grepl("#4D4D4D", def, fixed = TRUE))
    expect_true(grepl("#56B4E9", def, fixed = TRUE))

    session$setInputs(s_reference_colour = "#123456", s_comparison_colour = "#ABCDEF", generate = 2)
    expect_null(state$plot_error)
    alt <- svg_now()
    expect_true(grepl("#123456", alt, fixed = TRUE))
    expect_true(grepl("#ABCDEF", alt, fixed = TRUE))
    expect_false(grepl("#4D4D4D", alt, fixed = TRUE))
    expect_false(grepl("#56B4E9", alt, fixed = TRUE))
  })
})

test_that("uniprot_reviewed_only switches the Gene info UniProt link to Swiss-Prot only (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  # The panel also makes a live UniProt call (canonical-isoform check); fail it fast.
  # It shows the rendered figure's gene, so a (stubbed) figure is generated first.
  local_mocked_bindings(.uniprot_read_lines = function(url) stop("offline (test)"),
                        render_all = function(entries, region_params, settings, device_fn) invisible(NULL))
  cp <- effects_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    session$setInputs(view_mode = "compose", region_mode = "gene", gene_name = "Elf5", figure_tracks = "genes")
    session$setInputs(generate = 1)
    link <- function() {
      html <- as.character(output$gr_extras_ui$html)
      utils::URLdecode(sub(".*href=\"(https://www.uniprot.org/[^\"]*)\".*", "\\1", html))
    }
    def <- link()
    expect_match(def, "(gene:Elf5) AND (organism_id:10090)", fixed = TRUE)
    expect_false(grepl("reviewed", def, fixed = TRUE))
    session$setInputs(s_uniprot_reviewed_only = TRUE)
    expect_match(link(), "AND (reviewed:true)", fixed = TRUE)
  })
})

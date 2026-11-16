# The gene box's message (an error, or "Found") survives a trip to Coordinates and back:
# nothing re-runs the lookup on the way back, so a cleared one never returned.

test_that("Gene, Coordinates, then Gene again keeps the gene box's error or Found", {
  proj <- tempfile("pt-region-mode-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  suppressMessages(suppressWarnings(pt_project(proj, pt_example("mammary"))))
  json <- file.path(proj, "project_settings.json")
  app <- .pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))
  shiny::testServer(app, {
    box <- function() paste(as.character(output$gene_result_ui$html))
    session$setInputs(view_mode = "compose", region_mode = "gene", gene_name = "Elf5")
    expect_match(box(), "Found")
    session$setInputs(region_mode = "coords")
    session$setInputs(region_mode = "gene")
    expect_match(box(), "Found")
    # An invalid flank: the gene is cut away, an error instead of a region.
    session$setInputs(s_gene_flank_pct = -60)
    expect_match(box(), "leaves nothing of Elf5")
    session$setInputs(region_mode = "coords")
    session$setInputs(region_mode = "gene")
    expect_match(box(), "leaves nothing of Elf5")
    # A gene that isn't there, likewise.
    session$setInputs(s_gene_flank_pct = 15, gene_name = "Nofoundgene")
    expect_match(box(), "not found")
    session$setInputs(region_mode = "coords")
    session$setInputs(region_mode = "gene")
    expect_match(box(), "not found")
  })
})

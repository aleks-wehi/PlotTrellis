# The Gene / protein info panel (R/gene_info_panel.R) and the Figure settings layout,
# driving the real composer server headlessly (pattern of test-composer-settings.R).
# The panel sits on the Output page and follows the RENDERED figure's gene; its
# UniProt parts show a loading note only when their cache is cold. UniProt is never
# contacted and no disk cache is read or written: its cache, network reads and the
# protein-domain lookup are replaced in memory. render_all() is stubbed (the figure
# itself doesn't matter here).

# No project open (an earlier test file may have opened one): caches go beside the GTF.
.pt_set_cache_root(NULL)

gi_composer <- function() {
  proj <- tempfile("pt-gene-info-"); dir.create(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  list(proj = proj, app = suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))
}

# In-memory UniProt: the canonical-isoform cache, canned reads (counted: a search hit
# Q99999, then its 99 aa record, which matches no transcript), and a protein-domain
# lookup that is "cached" once it has run.
local_fake_uniprot <- function(env = parent.frame()) {
  mem <- new.env()
  mem$cache <- list(); mem$reads <- 0L; mem$pd_genes <- character(0); mem$pd_calls <- 0L
  local_mocked_bindings(
    render_all = function(entries, region_params, settings, device_fn) invisible(NULL),
    .uniprot_canonical_cache_load = function() mem$cache,
    .uniprot_canonical_cache_save = function(cache) mem$cache <- cache,
    .uniprot_read_lines = function(url) {
      mem$reads <- mem$reads + 1L
      if (grepl("/search?", url, fixed = TRUE)) "Q99999"
      else c("ID   TEST_MOUSE              Reviewed;          99 AA.",
             "SQ   SEQUENCE   99 AA;  11111 MW;  0000000000000000 CRC64;")
    },
    .pt_gi_pd_cached = function(gene, genome_label) gene %in% mem$pd_genes,
    gene_protein_domains = function(gene_symbol, ...) {
      mem$pd_calls <- mem$pd_calls + 1L
      mem$pd_genes <- union(mem$pd_genes, gene_symbol)
      list(status = "no_domains", gene = gene_symbol)
    },
    .env = env)
  mem
}

gene_figure <- function(gene) list(view_mode = "compose", region_mode = "gene", gene_name = gene,
                                   figure_tracks = "genes")
LOADING <- "Fetching UniProt annotations"

test_that("the panel follows the rendered figure's gene, not a later Compose change, and is absent for coordinates", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  local_fake_uniprot()
  cp <- gi_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  g <- GENOMES$mouse
  n_tx <- nrow(transcript_annotations("Nfatc1", TXDBS$mouse, g$orgdb, g$tx_name_lookup_path, g$tx_meta_lookup_path))
  shiny::testServer(cp$app, {
    html <- function() as.character(output$gr_extras_ui$html)
    do.call(session$setInputs, gene_figure("Nfatc1"))
    # A gene chosen on Compose but no figure yet: no panel.
    expect_identical(output$gene_info_gene, "")
    expect_null(output$gr_extras_ui)
    session$setInputs(generate = 1)
    expect_null(state$plot_error)
    expect_identical(output$gene_info_gene, "Nfatc1")
    # The transcript table: one <details> row per transcript of the gene.
    expect_equal(lengths(regmatches(html(), gregexpr("<details", html(), fixed = TRUE))), n_tx)
    expect_match(html(), ">Nfatc1-201<", fixed = TRUE)
    expect_gt(n_tx, 1)

    # Changing the gene on Compose doesn't change the panel.
    session$setInputs(gene_name = "Elf5")
    expect_identical(state$gene_symbol, "Elf5")
    expect_identical(output$gene_info_gene, "Nfatc1")
    expect_match(html(), "Nfatc1-201", fixed = TRUE)
    expect_false(grepl("Elf5", html(), fixed = TRUE))
    # Zoom re-renders the same figure's gene.
    session$setInputs(zoom_in_15x = 1)
    expect_identical(output$gene_info_gene, "Nfatc1")
    # A new figure for the new gene moves the panel.
    session$setInputs(generate = 2)
    expect_identical(output$gene_info_gene, "Elf5")
    expect_match(html(), "Elf5-201", fixed = TRUE)
    expect_false(grepl("Nfatc1", html(), fixed = TRUE))

    # A coordinates-only figure: no gene, no panel.
    session$setInputs(region_mode = "coords", region_chr = "chr2",
                      region_range = c(103300000L, 103500000L), generate = 3)
    expect_null(state$plot_error)
    expect_identical(output$gene_info_gene, "")
    expect_null(output$gr_extras_ui)
    expect_null(output$gr_protein_ui)
  })
})

test_that("UniProt parts show a loading note when their cache is cold, and fill in at once when warm", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  mem <- local_fake_uniprot()
  cp <- gi_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  # Cold: the first value sent shows the transcript table with a loading note; the
  # fetch runs after that flush, and the next flush fills the UniProt lines in.
  shiny::testServer(cp$app, {
    do.call(session$setInputs, gene_figure("Elf5"))
    session$setInputs(generate = 1)
    cold <- as.character(output$gr_extras_ui$html)
    expect_match(cold, LOADING, fixed = TRUE)
    expect_match(cold, "Elf5-201", fixed = TRUE)
    expect_match(as.character(output$gr_protein_ui$html), LOADING, fixed = TRUE)
    session$flushReact()
    filled <- as.character(output$gr_extras_ui$html)
    expect_false(grepl(LOADING, filled, fixed = TRUE))
    expect_false(grepl("Q99999", cold, fixed = TRUE))
    expect_match(filled, "canonical is Q99999", fixed = TRUE)
    expect_match(as.character(output$gr_protein_ui$html), "No annotated domain/region/motif features", fixed = TRUE)
    expect_false(grepl(LOADING, as.character(output$gr_protein_ui$html), fixed = TRUE))
  })
  expect_equal(mem$reads, 2L)  # one search + one record, then cached
  expect_true("org.Mm.eg.db|Elf5" %in% names(mem$cache))
  expect_equal(mem$pd_calls, 1L)

  # Warm (a new session, so nothing is held in memory): filled from the first value on.
  shiny::testServer(cp$app, {
    do.call(session$setInputs, gene_figure("Elf5"))
    session$setInputs(generate = 1)
    warm <- as.character(output$gr_extras_ui$html)
    expect_false(grepl(LOADING, warm, fixed = TRUE))
    expect_match(warm, "canonical is Q99999", fixed = TRUE)
    expect_match(as.character(output$gr_protein_ui$html), "No annotated domain/region/motif features", fixed = TRUE)
  })
  expect_equal(mem$reads, 2L)  # no new network read
})

test_that("a transcript_annotations() failure shows an error line instead of an empty table", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  local_fake_uniprot()
  local_mocked_bindings(transcript_annotations = function(...) stop("lookup file missing (test)"))
  cp <- gi_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    do.call(session$setInputs, gene_figure("Elf5"))
    session$setInputs(generate = 1)
    html <- as.character(output$gr_extras_ui$html)
    expect_match(html, "load the transcript annotations:", fixed = TRUE)
    expect_match(html, "lookup file missing (test)", fixed = TRUE)
    expect_false(grepl("<details", html, fixed = TRUE))
    expect_match(html, "UniProt: Elf5", fixed = TRUE)  # the link-out still works
  })
})

test_that("the panel is on the Output page and the title field is first in Figure settings", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- gi_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  html <- as.character(get("ui", envir = environment(cp$app$serverFuncSource())))
  at <- function(s) regexpr(s, html, fixed = TRUE)[[1]]
  # Gene / protein info: inside the Output section, after the downloads, before the R
  # script, shown only for a gene figure.
  expect_gt(at('id="gene_protein_info-body"'), at('id="output-section"'))
  expect_gt(at('id="gene_protein_info-body"'), at('id="dl_rcode"'))
  expect_lt(at('id="gene_protein_info-body"'), at('id="rcode-body"'))
  expect_match(html, 'data-display-if="output.gene_info_gene"', fixed = TRUE)
  # Until its first value arrives, the table area says it is loading.
  expect_match(html, '<div id="gr_extras_ui" class="shiny-html-output">\\s*<div[^>]*>Loading transcript annotations...</div>')
  # Title: inside the Figure settings panel, above the grouped settings, before Generate.
  expect_gt(at('id="plot_title"'), at('id="figure_settings-body"'))
  expect_lt(at('id="plot_title"'), at('id="display_settings_ui"'))
  expect_lt(at('id="plot_title"'), at('id="save_settings"'))
})

test_that("every settings group is a collapsed subsection and every s_<key> input is still rendered", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- gi_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    html <- output$display_settings_ui$html
    groups <- regmatches(html, gregexpr('<details class="pt-settings-group"[^>]*>\\s*<summary[^>]*>[^<]*', html))[[1]]
    expect_equal(sub(".*>", "", groups),
                 c("Layout &amp; display", "Coverage tracks (BigWig: ATAC / ChIP-seq / CUT&amp;RUN / RNA-seq)",
                   "Hi-C", "Junction / arc tracks", "Gene view", "Comparison (diff) tracks"))
    expect_false(any(grepl(" open", groups, fixed = TRUE)))  # all collapsed
    for (k in setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
      expect_true(grepl(sprintf('id="s_%s"', k), html, fixed = TRUE), info = k)
    # The title is not a setting: setting it leaves the settings clean.
    session$setInputs(plot_title = "My title")
    expect_false(settings_dirty())
  })
})

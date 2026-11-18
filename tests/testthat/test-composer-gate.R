# The composer's stage gating (R/app-gate.R) on its own: the region and conditions
# blocks are replaced by reactive values the test sets, the track items are hand-made
# (helper-composer-fakes.R), and the lock messages are recorded by the mock session.

gate_app <- function(ctx = fake_ctx()) function(input, output, session) {
  state <- reactiveValues(gene_region = NULL, gene_symbol = NULL)
  # What region_valid() gives; NULL makes it stop silently, as before any region is set.
  rv <- reactiveVal(NULL)
  ready <- reactiveVal(FALSE)
  gt <- .pt_generate_gate_server(input, output, session, state, ctx,
                                 region_valid = reactive({ req(rv()); rv() }),
                                 has_conditions_ready = reactive(ready()))
}

VALID <- list(valid = TRUE, chrom = "chr1", start = 100L, end = 2000L)
GENE_RGN <- list(chrom = "chr4", start = 44700000L, end = 44900000L)

# The `locked` values of one lock message type, in the order sent.
locks <- function(rec, type)
  vapply(Filter(function(m) identical(m$type, type), rec$custom), function(m) m$message$locked, logical(1))

generate_enabled <- function(html) !grepl("disabled", html, fixed = TRUE) && grepl('id="generate"', html, fixed = TRUE)

test_that("Coordinates resolves to the valid region with no gene name", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    session$setInputs(region_mode = "coords")
    expect_null(gt$resolved_region())
    rv(VALID)
    expect_identical(gt$resolved_region(),
                     list(rgn = list(chrom = "chr1", start = 100L, end = 2000L), gene_name = NA_character_))
    rv(list(valid = FALSE, msg = "End must be > start"))
    expect_null(gt$resolved_region())
  }, session = rec$session)
})

test_that("Coordinates ignores a found gene, and a gene region ignores the coordinates", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    state$gene_region <- GENE_RGN
    state$gene_symbol <- "Pax5"
    session$setInputs(region_mode = "coords")
    expect_null(gt$resolved_region())
    rv(VALID)
    session$setInputs(region_mode = "gene")
    expect_identical(gt$resolved_region(), list(rgn = GENE_RGN, gene_name = "Pax5"))
  }, session = rec$session)
})

test_that("a gene region is named by the resolved symbol, else the gene box", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    session$setInputs(region_mode = "gene", gene_name = "ENSMUSG00000014030")
    expect_null(gt$resolved_region())
    state$gene_region <- GENE_RGN
    expect_identical(gt$resolved_region()$gene_name, "ENSMUSG00000014030")
    state$gene_symbol <- "Pax5"
    expect_identical(gt$resolved_region()$gene_name, "Pax5")
  }, session = rec$session)
})

test_that("with no region mode chosen, a gene region still resolves", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    state$gene_region <- GENE_RGN
    # Any mode other than "coords" (here none yet) reads the gene lookup's region.
    expect_identical(gt$resolved_region()$rgn, GENE_RGN)
    expect_identical(gt$resolved_region()$gene_name, NULL)
  }, session = rec$session)
})

test_that("has_mode is whether a sample-count mode is chosen", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    expect_false(gt$has_mode())
    session$setInputs(view_mode = "single")
    expect_true(gt$has_mode())
  }, session = rec$session)
})

test_that("the condition lock opens once a mode and its conditions are ready", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    session$flushReact()
    expect_identical(locks(rec, "condition_lock"), TRUE)
    session$setInputs(view_mode = "two")
    ready(TRUE); session$flushReact()
    expect_identical(tail(locks(rec, "condition_lock"), 1), FALSE)
    ready(FALSE); session$flushReact()
    expect_identical(tail(locks(rec, "condition_lock"), 1), TRUE)
  }, session = rec$session)
})

test_that("the region lock also waits for a resolved region", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    session$setInputs(view_mode = "compose", region_mode = "coords")
    ready(TRUE); session$flushReact()
    expect_true(all(locks(rec, "region_lock")))
    expect_identical(tail(locks(rec, "condition_lock"), 1), FALSE)
    rv(VALID); session$flushReact()
    expect_identical(tail(locks(rec, "region_lock"), 1), FALSE)
    # Switching to Gene name with no gene found locks the picker again.
    session$setInputs(region_mode = "gene")
    expect_identical(tail(locks(rec, "region_lock"), 1), TRUE)
  }, session = rec$session)
})

test_that("invalid_figure_tracks lists broken tracks, and gene-only tracks in Coordinates", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    expect_identical(gt$invalid_figure_tracks(), character(0))
    session$setInputs(figure_tracks = c("genes", "bw_atac_A", "bw_atac_B", "tpm_usage_Quant_A", "transcripts"))
    expect_identical(gt$invalid_figure_tracks(), "bw_atac_B")
    session$setInputs(region_mode = "coords")
    expect_identical(gt$invalid_figure_tracks(), c("bw_atac_B", "tpm_usage_Quant_A", "transcripts"))
    # A key with no item is skipped, not reported.
    session$setInputs(figure_tracks = c("no_such_key", "genes"))
    expect_identical(gt$invalid_figure_tracks(), character(0))
  }, session = rec$session)
})

test_that("Generate is enabled only with a mode, conditions, a region, tracks and no invalid track", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    expect_false(generate_enabled(ui_html(output$generate_ui)))
    session$setInputs(view_mode = "single", region_mode = "coords", figure_tracks = "genes")
    rv(VALID)
    expect_false(generate_enabled(ui_html(output$generate_ui)))   # conditions not ready
    ready(TRUE); session$flushReact()
    expect_true(generate_enabled(ui_html(output$generate_ui)))
    session$setInputs(figure_tracks = character(0))
    expect_false(generate_enabled(ui_html(output$generate_ui)))   # no tracks
    session$setInputs(figure_tracks = c("genes", "transcripts"))
    expect_false(generate_enabled(ui_html(output$generate_ui)))   # a gene-only track in Coordinates
    rv(list(valid = FALSE, msg = "x")); session$setInputs(figure_tracks = "genes")
    expect_false(generate_enabled(ui_html(output$generate_ui)))   # no region
  }, session = rec$session)
})

test_that("Generate stays disabled with no mode even when everything else is ready", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    session$setInputs(region_mode = "coords", figure_tracks = "genes")
    rv(VALID); ready(TRUE); session$flushReact()
    expect_false(generate_enabled(ui_html(output$generate_ui)))
    session$setInputs(view_mode = "compose")
    expect_true(generate_enabled(ui_html(output$generate_ui)))
  }, session = rec$session)
})

test_that("the warning beside Generate names the invalid tracks by label", {
  rec <- recording_session()
  shiny::testServer(gate_app(), {
    expect_null(output$generate_warning_ui)
    session$setInputs(region_mode = "coords", figure_tracks = c("bw_atac_B", "genes", "transcripts"))
    html <- ui_html(output$generate_warning_ui)
    expect_match(html, '<div class="alert alert-warning" style="margin-bottom:0;">', fixed = TRUE)
    expect_match(html, "<strong>Remove or fix before generating: </strong>", fixed = TRUE)
    expect_match(html, "atac \u00B7 B, Transcripts (filtered)", fixed = TRUE)
  }, session = rec$session)
})

test_that("an invalid track with no label is named by its key", {
  ctx <- fake_ctx()
  ctx$TRACK_ITEMS$bw_atac_B$label <- NULL
  rec <- recording_session()
  shiny::testServer(gate_app(ctx), {
    session$setInputs(figure_tracks = "bw_atac_B")
    expect_match(ui_html(output$generate_warning_ui), "</strong>\\s*bw_atac_B\\s*</div>")
  }, session = rec$session)
})

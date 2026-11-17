# The composer's gene lookup (R/composer-gene.R) on its own: .pt_resolve_gene_region() is
# replaced by a lookup table, so no genome is needed. Each lookup's arguments and the gene
# box updates (the canonical symbol echoed back) are recorded.

gene_app <- function(ctx = fake_ctx(), saved = list(gene_flank_pct = 15)) function(input, output, session) {
  state <- reactiveValues(gene_region = NULL, gene_error = NULL, gene_note = NULL, gene_symbol = NULL)
  saved_cfg <- reactiveVal(saved)
  .pt_gene_lookup_server(input, output, session, state, ctx, saved_cfg)
}

PAX5_RGN <- list(chrom = "chr4", start = 44700000L, end = 44900000L)
ENS_RGN  <- list(chrom = "chr4", start = 44710000L, end = 44890000L)

# A stand-in for .pt_resolve_gene_region(): Pax5 in any case, its Ensembl id (with its own
# region and a note), and a flank-error case; anything else is not found.
fake_resolve <- function(query, txdb, flank_pct, orgdb_pkg, all_genes, genome_label) {
  if (tolower(query) == "pax5") list(region = PAX5_RGN, symbol = "Pax5", error = NULL)
  else if (query == "ENSMUSG00000014030")
    list(region = ENS_RGN, symbol = "Pax5", error = NULL, note = "ENSMUSG00000014030 is also named Pax5b")
  else if (query == "elf5") list(region = NULL, symbol = "Elf5", error = "A gene flank of -60% leaves nothing of Elf5 to show: use a larger flank.")
  else list(region = NULL, symbol = NA_character_,
            error = paste0("'", query, "' not found in ", genome_label, ": enter a gene symbol (e.g. Pax5) or an Ensembl gene ID"))
}

# Runs `code` against the gene lookup server (it sees session, output and state); returns
# the lookups made (their arguments) and the gene box updates.
drive_gene <- function(code, ctx = fake_ctx(), saved = list(gene_flank_pct = 15)) {
  code <- rlang::enquo(code)
  rec <- new.env(); rec$lookups <- list(); rec$updates <- list()
  local_mocked_bindings(.pt_resolve_gene_region = function(query, txdb, flank_pct, orgdb_pkg, all_genes, genome_label) {
    rec$lookups[[length(rec$lookups) + 1]] <- list(query = query, txdb = txdb, flank_pct = flank_pct, orgdb_pkg = orgdb_pkg,
                                                   all_genes = all_genes, genome_label = genome_label)
    fake_resolve(query, txdb, flank_pct, orgdb_pkg, all_genes, genome_label)
  })
  local_mocked_bindings(updateSelectizeInput = function(session, inputId, label = NULL, choices = NULL,
                                                        selected = NULL, options = list(), server = FALSE)
    rec$updates[[length(rec$updates) + 1]] <- list(id = inputId, choices = choices, selected = selected, server = server))
  rlang::inject(shiny::testServer(gene_app(ctx, saved), !!code))
  rec
}

queries <- function(rec) vapply(rec$lookups, `[[`, "", "query")

test_that("an empty gene box looks nothing up and shows nothing", {
  r <- drive_gene({
    session$flushReact()
    expect_null(state$gene_region)
    expect_null(output$gene_result_ui)
    session$setInputs(gene_name = "   ")
    expect_null(state$gene_region)
  })
  expect_length(r$lookups, 0)
})

test_that("a found gene sets the region and symbol, passing the context's genome to the lookup", {
  r <- drive_gene({
    session$setInputs(gene_name = "Pax5")
    expect_identical(state$gene_region, PAX5_RGN)
    expect_identical(state$gene_symbol, "Pax5")
    expect_null(state$gene_error)
    expect_null(state$gene_note)
  })
  expect_identical(r$lookups[[1]], list(query = "Pax5", txdb = "fake-txdb", flank_pct = 0.15, orgdb_pkg = "org.Fake.eg.db",
                                        all_genes = c("Elf5", "Pax5"), genome_label = "mm10"))
  # The symbol typed is already canonical: no echo.
  expect_length(r$updates, 0)
})

test_that("the result line gives the region, its length and the note", {
  drive_gene({
    session$setInputs(gene_name = "Pax5")
    expect_identical(ui_html(output$gene_result_ui),
                     '<div class="alert alert-success">\n  <strong>Found: </strong>\n  chr4:44.7Mb \u2013 44.9Mb (200Kb)\n</div>')
    session$setInputs(gene_name = "ENSMUSG00000014030")
    expect_match(ui_html(output$gene_result_ui), "chr4:44.71Mb \u2013 44.89Mb (180Kb)\n  <br/>\n  ENSMUSG00000014030 is also named Pax5b", fixed = TRUE)
  })
})

test_that("a gene not found clears the region and shows the lookup's error", {
  r <- drive_gene({
    session$setInputs(gene_name = "Pax5")
    session$setInputs(gene_name = "Nope1")
    expect_null(state$gene_region)
    expect_identical(state$gene_symbol, NA_character_)
    expect_identical(state$gene_error, "'Nope1' not found in mm10: enter a gene symbol (e.g. Pax5) or an Ensembl gene ID")
    expect_identical(ui_html(output$gene_result_ui),
                     "<div class=\"alert alert-danger\">'Nope1' not found in mm10: enter a gene symbol (e.g. Pax5) or an Ensembl gene ID</div>")
  })
  expect_length(r$updates, 0)
})

test_that("a wrong-case symbol is echoed back once as the canonical symbol, with no second lookup", {
  r <- drive_gene({
    session$setInputs(gene_name = "pax5")
    expect_identical(state$gene_symbol, "Pax5")
    # The browser then reports the echoed value: skipped, the region stays.
    session$setInputs(gene_name = "Pax5")
    expect_identical(state$gene_region, PAX5_RGN)
  })
  expect_identical(queries(r), "pax5")
  expect_identical(r$updates, list(list(id = "gene_name", choices = c("Elf5", "Pax5"), selected = "Pax5", server = TRUE)))
})

test_that("an Ensembl id keeps its own region after the echo and after a flank change", {
  r <- drive_gene({
    session$setInputs(gene_name = "ENSMUSG00000014030")
    expect_identical(state$gene_region, ENS_RGN)
    session$setInputs(gene_name = "Pax5")      # the echo
    expect_identical(state$gene_region, ENS_RGN)
    session$setInputs(s_gene_flank_pct = 30)   # re-resolves the id typed, not the symbol
    expect_identical(state$gene_region, ENS_RGN)
  })
  expect_identical(queries(r), c("ENSMUSG00000014030", "ENSMUSG00000014030"))
  expect_identical(vapply(r$lookups, `[[`, 0, "flank_pct"), c(0.15, 0.30))
})

test_that("clearing the box clears the whole lookup and the remembered query", {
  r <- drive_gene({
    session$setInputs(gene_name = "ENSMUSG00000014030")
    session$setInputs(gene_name = "Pax5")
    session$setInputs(gene_name = "")
    expect_null(state$gene_region); expect_null(state$gene_symbol)
    expect_null(state$gene_error); expect_null(state$gene_note)
    expect_null(output$gene_result_ui)
    # Typed afresh, the symbol is looked up as itself.
    session$setInputs(gene_name = "Pax5")
    expect_identical(state$gene_region, PAX5_RGN)
  })
  expect_identical(queries(r), c("ENSMUSG00000014030", "Pax5"))
})

test_that("the box is trimmed before the lookup", {
  r <- drive_gene(session$setInputs(gene_name = "  Pax5 "))
  expect_identical(queries(r), "Pax5")
  # The trimmed text is the canonical symbol, so there is no echo.
  expect_length(r$updates, 0)
})

test_that("a symbol found with no region (a flank error) is not echoed", {
  r <- drive_gene({
    session$setInputs(gene_name = "elf5")
    expect_null(state$gene_region)
    expect_identical(state$gene_symbol, "Elf5")
    expect_match(state$gene_error, "leaves nothing of Elf5", fixed = TRUE)
  })
  expect_length(r$updates, 0)
})

test_that("the flank comes from the panel's box, else the saved settings, else the default", {
  flank <- function(inputs, saved) {
    r <- drive_gene(do.call(session$setInputs, c(list(gene_name = "Pax5"), inputs)), saved = saved)
    r$lookups[[length(r$lookups)]]$flank_pct
  }
  expect_identical(flank(list(s_gene_flank_pct = 40), list(gene_flank_pct = 15)), 0.40)
  expect_identical(flank(list(), list(gene_flank_pct = 25)), 0.25)
  expect_identical(flank(list(s_gene_flank_pct = "abc"), list(gene_flank_pct = 25)), 0.25)
  expect_identical(flank(list(s_gene_flank_pct = NA), list(gene_flank_pct = "25")), 0.25)
  expect_identical(flank(list(s_gene_flank_pct = c(10, 20)), list(gene_flank_pct = 25)), 0.25)
  expect_identical(flank(list(s_gene_flank_pct = Inf), list(gene_flank_pct = "wide")), 0.15)
  expect_identical(flank(list(), list()), 0.15)
  # 0 and negative values are passed on; the lookup decides what they leave.
  expect_identical(flank(list(s_gene_flank_pct = 0), list(gene_flank_pct = 15)), 0)
  expect_identical(flank(list(s_gene_flank_pct = -50), list(gene_flank_pct = 15)), -0.5)
})

test_that("an error result takes the place of an earlier region in the result line", {
  drive_gene({
    session$flushReact()  # the lookup's first run, with an empty box
    state$gene_region <- PAX5_RGN
    state$gene_error <- "Something failed"
    session$flushReact()
    expect_identical(ui_html(output$gene_result_ui), '<div class="alert alert-danger">Something failed</div>')
  })
})

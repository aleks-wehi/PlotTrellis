# The composer's gene lookup (R/app-main.R): resolves the gene box to a region in
# `state`, echoes the canonical symbol back into the box, and shows the result
# line. The flank comes from the Figure settings input, else the saved settings.
.pt_gene_lookup_server <- function(input, output, session, state, ctx,
                                   saved_cfg) {
  # ── Region: gene lookup (shared by all 3 modes) ─────────────────────────────
  # The symbol this observer itself just put in the box, so the re-run that
  # update triggers keeps a typed Ensembl ID's own region rather than
  # re-resolving the symbol to whichever of its IDs comes first.
  gene_echo <- NULL
  # The text the current gene was resolved from, so a flank-only re-run
  # (Figure settings) re-resolves a typed Ensembl ID, not its echoed symbol.
  gene_query <- NULL
  observe({
    sym <- trimws(input$gene_name %||% "")
    # Read first, so the echo run below still depends on it. A non-numeric box
    # or saved value (e.g. a hand-edited project file) falls back to the
    # default, rather than an NA flank erroring here and ending the session.
    pct <- Find(function(x) length(x) == 1 && is.finite(x),
                lapply(list(input$s_gene_flank_pct, saved_cfg()$gene_flank_pct),
                       function(x) suppressWarnings(as.numeric(x))),
                nomatch = .pt_default_settings()$gene_flank_pct) / 100
    if (identical(sym, gene_echo)) { gene_echo <<- NULL; return() }
    gene_echo <<- NULL
    if (nchar(sym) == 0) {
      # An empty box is no gene: clear the whole lookup, so Generate has no
      # region to draw.
      state$gene_symbol <- NULL; state$gene_region <- NULL
      state$gene_error <- NULL; state$gene_note <- NULL
      gene_query <<- NULL
      return()
    }
    q <- if (!is.null(gene_query) &&
             identical(sym, isolate(state$gene_symbol)))
      gene_query else sym
    gene_query <<- q
    res <- .pt_resolve_gene_region(q, ctx$txdb, pct, ctx$ORGDB_PKG,
                                   ctx$ALL_GENES, ctx$GENOME$genome_label)
    state$gene_region <- res$region
    state$gene_error  <- res$error
    state$gene_symbol <- res$symbol
    state$gene_note   <- res$note
    # An Ensembl ID or wrong-case symbol resolved: show the canonical symbol,
    # which is what every track and panel below looks genes up by.
    if (!is.null(res$region) && !identical(res$symbol, sym)) {
      gene_echo <<- res$symbol
      # A gene named only in the GTF is not among the OrgDb's symbols.
      updateSelectizeInput(session, "gene_name",
                           choices = union(ctx$ALL_GENES, res$symbol),
                           server = TRUE, selected = res$symbol)
    }
  })

  output$gene_result_ui <- renderUI({
    if (!is.null(state$gene_error))
      div(class = "alert alert-danger", state$gene_error)
    else if (!is.null(state$gene_region)) {
      r <- state$gene_region
      div(class = "alert alert-success",
          strong("Found: "),
          paste0(r$chrom, ":", .pt_format_coord(r$start), " \u2013 ",
                 .pt_format_coord(r$end),
                 " (", .pt_format_bp(r$end - r$start), ")"),
          if (!is.null(state$gene_note)) tagList(br(), state$gene_note))
    }
  })

  invisible(NULL)
}

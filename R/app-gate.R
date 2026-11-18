# The composer's stage gating (R/app-main.R): the resolved region (gene or
# coordinates), whether a mode is chosen, the condition and region locks sent to
# the page, the Figure tracks that cannot draw, and the Generate button with its
# warning. `region_valid` and `has_conditions_ready` are the region and
# conditions blocks' reactives.
.pt_generate_gate_server <- function(input, output, session, state, ctx,
                                     region_valid, has_conditions_ready) {
  # ── Region resolution: identical across all 3 sample-count modes ────────────
  # Coordinates is selectable in every mode, so which region resolves depends
  # only on region_mode -- never on view_mode. gene_name is only meaningful in
  # 'gene' region_mode; tracks that need a gene symbol (GENE_REQUIRED_TYPES)
  # can't draw from coordinates in any mode, which the picker's red-boxing
  # handles, not this resolution logic.
  resolved_region <- reactive({
    if (!is.null(input$region_mode) && input$region_mode == "coords") {
      v <- tryCatch(region_valid(), error = function(e) NULL)
      if (!is.null(v) && isTRUE(v$valid))
        list(rgn = list(chrom = v$chrom, start = v$start, end = v$end),
             gene_name = NA_character_)
      else NULL
    } else if (!is.null(state$gene_region)) {
      list(rgn = state$gene_region,
           gene_name = state$gene_symbol %||% input$gene_name)
    } else NULL
  })

  # Whether a sample-count mode has been chosen at all (view_mode starts
  # unselected) -- Generate and the locks below all stay gated on this in
  # addition to conditions/region, so the intended stage order (mode ->
  # conditions -> region/gene -> tracks -> generate) holds whatever else is
  # already set.
  has_mode <- reactive({ !is.null(input$view_mode) })

  # region_mode/gene_name lock until mode AND conditions are both resolved,
  # so a gene can't be searched before Reference/Comparison (or Condition)
  # are chosen. Mirrors the 'region_lock' message
  # below (same server-driven pattern), consumed by setConditionGateLocked().
  observe({
    session$sendCustomMessage("condition_lock",
                              list(locked = !(has_mode() &&
                                              has_conditions_ready())))
  })

  # The track picker (both Available and Figure tracks) locks whenever
  # mode+conditions+region aren't all ready -- covers the "nothing chosen yet"
  # case and the mid-composition one (switching to Coordinates with no
  # coordinates entered locks the existing Figure tracks too, since it's driven
  # by the same single flag rather than two separate lock layers). The
  # track-pool populator (R/app-conditions.R) also waits for the
  # conditions.
  observe({
    session$sendCustomMessage("region_lock",
                              list(locked = !(has_mode() &&
                                              has_conditions_ready() &&
                                              !is.null(resolved_region()))))
  })

  # Figure-tracks keys that can't draw as things stand -- same two conditions
  # the client already red-boxes (the track-gene-required/track-broken CSS
  # classes, applyGeneRequiredFlag()/applyBrokenFlag() in inst/www/composer.js),
  # reimplemented server-side so Generate is disabled while one is there.
  invalid_figure_tracks <- reactive({
    keys <- input$figure_tracks %||% character(0)
    if (length(keys) == 0) return(character(0))
    is_coords <- identical(input$region_mode, "coords")
    Filter(function(k) {
      item <- ctx$TRACK_ITEMS[[k]]
      !is.null(item) && (isTRUE(item$broken) ||
                         (is_coords && item$type %in% ctx$GENE_REQUIRED_TYPES))
    }, keys)
  })

  # ── Generate button (disabled until mode + conditions + region + tracks,
  # and no invalid track sitting in Figure tracks) ───────────────────────────
  output$generate_ui <- renderUI({
    labels     <- input$figure_tracks %||% character(0)
    has_tracks <- length(labels) > 0
    has_region <- !is.null(resolved_region())
    has_valid_tracks <- length(invalid_figure_tracks()) == 0
    .pt_generate_button_ui("generate",
                           has_mode() && has_conditions_ready() &&
                             has_region && has_tracks && has_valid_tracks)
  })

  # Message beside Generate naming the red-boxed track(s) blocking it -- mirrors
  # the .track-gene-required/.track-broken title tooltip text
  # (applyGeneRequiredFlag()/applyBrokenFlag() in inst/www/composer.js) but as
  # visible text, since a disabled button with no explanation is a cryptic dead
  # end.
  output$generate_warning_ui <- renderUI({
    bad <- invalid_figure_tracks()
    if (length(bad) == 0) return(NULL)
    names_str <- paste(vapply(bad,
                              function(k) ctx$TRACK_ITEMS[[k]]$label %||% k,
                              character(1)),
                       collapse = ", ")
    div(class = "alert alert-warning", style = "margin-bottom:0;",
        strong("Remove or fix before generating: "), names_str)
  })

  list(resolved_region = resolved_region, has_mode = has_mode,
       invalid_figure_tracks = invalid_figure_tracks)
}

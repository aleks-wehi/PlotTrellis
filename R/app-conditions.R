# The composer's condition selection (R/app-main.R): Reference/Comparison mutual
# exclusion, which conditions the current mode has picked, and the track pool
# populator that re-scopes the picker when the mode or conditions change.
# Returns the reactives the gating reads.
.pt_conditions_server <- function(input, output, session, ctx) {
  # ── Two-sample reference/comparison mutual exclusion ──
  # Both boxes always offer every condition; the only guard is against picking
  # the same one in both -- sel_prev tracks each box's last-settled value so we
  # can tell which box just changed (vs. one only just relabeled by this same
  # observer) and clear the OTHER box on a collision. ignoreInit = TRUE since
  # the initial choices/markers are already baked into the UI; both start ""
  # (neither box preselects a condition), matching
  # sel_reference/sel_comparison's own defaults.
  sel_prev <- reactiveValues(control = "", treatment = "")

  observeEvent(list(input$sel_reference, input$sel_comparison), {
    ctrl <- input$sel_reference    %||% ""
    trt  <- input$sel_comparison  %||% ""

    if (!identical(ctrl, sel_prev$control) && nchar(ctrl) > 0 &&
        identical(ctrl, trt)) {
      # control just changed to collide with treatment -- clear treatment
      trt <- ""
    } else if (!identical(trt, sel_prev$treatment) && nchar(trt) > 0 &&
               identical(trt, ctrl)) {
      # treatment just changed to collide with control -- clear control
      ctrl <- ""
    }

    sel_prev$control   <- ctrl
    sel_prev$treatment <- trt

    updateSelectizeInput(session, "sel_reference",
      choices = .pt_condition_choices(ctx$AVAILABLE_CONDITIONS),
      selected = ctrl)
    updateSelectizeInput(session, "sel_comparison",
      choices = .pt_condition_choices(ctx$AVAILABLE_CONDITIONS), selected = trt)
  }, ignoreInit = TRUE)

  # ── Which condition(s) are currently selected, mode-dependent ──
  selected_conditions <- reactive({
    if (identical(input$view_mode, "two"))
      c(input$sel_reference, input$sel_comparison)
    else input$sel_condition
  })

  # Whether the condition-selection stage is fully resolved for the current
  # mode: Compose has no condition dropdown at all (the track picker is its own
  # sample selector), so it's trivially ready; single/two need their dropdown(s)
  # actually filled in, since none of sel_condition/
  # sel_reference/sel_comparison preselect a condition. Shared by the
  # region/gene lock, the track-pool populator, and Generate (R/app-gate.R)
  # so all three agree on the same stage boundary.
  has_conditions_ready <- reactive({
    if (identical(input$view_mode, "compose")) TRUE
    else if (identical(input$view_mode, "two"))
      nchar(input$sel_reference %||% "") > 0 &&
        nchar(input$sel_comparison %||% "") > 0
    else if (identical(input$view_mode, "single"))
      nchar(input$sel_condition %||% "") > 0
    else FALSE
  })

  # Re-scopes and re-defaults the track picker when the mode or the selected
  # condition(s) change, never on a region/gene change (that would wipe what the
  # user has placed; the picker is locked by region_lock instead).
  # figure_pool_key is the (mode, conditions) the pool was last built from: a
  # trigger that leaves it unchanged is skipped. Runs at startup too (no
  # ignoreInit): the picker's static markup holds Multi's pool, which the other
  # modes must replace.
  figure_pool_key <- reactiveVal(NULL)

  observeEvent(list(input$view_mode, input$sel_condition,
                    input$sel_reference, input$sel_comparison), {
    if (identical(input$view_mode, "compose")) {
      key   <- "compose"
      split <- list(
        available = .pt_track_labels(ctx$COMPOSE_AVAILABLE_KEYS,
                                     ctx$TRACK_ITEMS),
        figure    = list(
          genes = HTML(.pt_format_track_label(ctx$TRACK_ITEMS[["genes"]]$label))
        )
      )
      if (identical(ctx$START_VIEW$view_mode, "compose"))
        split <- list(
          available = .pt_track_labels(
            setdiff(ctx$COMPOSE_AVAILABLE_KEYS,
                    setdiff(ctx$START_VIEW$tracks, "genes")),
            ctx$TRACK_ITEMS
          ),
          figure    = .pt_track_labels(ctx$START_VIEW$tracks, ctx$TRACK_ITEMS)
        )
    } else {
      conds <- selected_conditions()
      if (length(conds) == 0 || anyNA(conds) || any(nchar(conds) == 0)) {
        key   <- "not_ready"
        split <- list(available = list(), figure = list())
      } else {
        key   <- paste(input$view_mode, paste(conds, collapse = "|"),
                       sep = "::")
        split <- .pt_track_split(conds, ctx$TRACK_ITEMS, ctx$START_VIEW)
      }
    }
    if (identical(key, figure_pool_key())) return()
    figure_pool_key(key)
    update_rank_list("available_tracks", labels = split$available,
                     session = session)
    update_rank_list("figure_tracks",    labels = split$figure,
                     session = session)
  })

  # Until the sample(s) are picked both track lists are empty: the Available
  # list says what to pick first. The text goes in a CSS variable, set only
  # before the pick, that plottrellis.css shows in an empty Available list.
  output$tracks_hint <- renderUI({
    if (has_conditions_ready()) return(NULL)
    msg <- switch(
      input$view_mode %||% "",
      single = "Pick a condition above to see its tracks.",
      two    = "Pick a reference and a comparison above to see their tracks.",
      "Pick the number of samples above to start."
    )
    tags$style(sprintf('#tracks-block { --pt-tracks-hint: "%s"; }', msg))
  })

  list(sel_prev = sel_prev, selected_conditions = selected_conditions,
       has_conditions_ready = has_conditions_ready,
       figure_pool_key = figure_pool_key)
}

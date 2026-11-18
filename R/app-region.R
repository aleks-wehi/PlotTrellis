# The composer's coordinates block (R/app-main.R): chromosome select, range slider
# and the start/end text boxes kept in step with it, their parse errors in
# `state`, and the region validity the gating reads. `chr_lengths` is the
# genome's named chromosome length vector.
.pt_region_server <- function(input, output, session, state, chr_lengths) {
  # ── Region: chromosome selector ─────────────────────────────────────────────
  output$chr_select_ui <- renderUI({
    chrs <- .pt_sort_chrs(names(chr_lengths))
    selectizeInput("region_chr", "Chromosome:", choices = chrs,
                   selected = chrs[1])
  })

  chr_len_r <- reactive({
    req(input$region_chr)
    .pt_chr_length(chr_lengths, input$region_chr, 200000000L)
  })

  output$range_slider_ui <- renderUI({
    len <- chr_len_r()
    sliderInput("region_range", "Region:",
                min = 1L, max = len, value = c(1L, len),
                step = 1000L, width = "100%", sep = ",")
  })

  output$start_num_ui <- renderUI({
    textInput("region_start_text", "Start:",
              value = .pt_format_coord(1L), width = "100%",
              placeholder = "e.g. 17Mb, 15.5kb, 17000000")
  })

  output$end_num_ui <- renderUI({
    textInput("region_end_text", "End:",
              value = .pt_format_coord(chr_len_r()), width = "100%",
              placeholder = "e.g. 17.49Mb, 17490kb")
  })

  # Slider → text: update display when slider moves
  observeEvent(input$region_range, {
    req(input$region_range)
    updateTextInput(session, "region_start_text",
                    value = .pt_format_coord(input$region_range[1]))
    updateTextInput(session, "region_end_text",
                    value = .pt_format_coord(input$region_range[2]))
  }, ignoreInit = TRUE)

  # Text → slider: parse on blur/Enter (JS suppresses auto-fire; only explicit
  # events reach here)
  observeEvent(input$region_start_text, {
    txt <- input$region_start_text
    if (is.null(txt) || nchar(trimws(txt)) == 0) return()
    req(input$region_range)
    v <- .pt_parse_bp(txt)
    if (is.na(v)) {
      state$start_parse_err <-
        "Cannot parse start. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)"
      return()
    }
    state$start_parse_err <- NULL
    if (!isTRUE(input$region_range[1] == v))
      updateSliderInput(session, "region_range",
                        value = c(v, input$region_range[2]))
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$region_end_text, {
    txt <- input$region_end_text
    if (is.null(txt) || nchar(trimws(txt)) == 0) return()
    req(input$region_range)
    v <- .pt_parse_bp(txt)
    if (is.na(v)) {
      state$end_parse_err <-
        "Cannot parse end. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)"
      return()
    }
    state$end_parse_err <- NULL
    if (!isTRUE(input$region_range[2] == v))
      updateSliderInput(session, "region_range",
                        value = c(input$region_range[1], v))
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  region_valid <- reactive({
    req(input$region_chr, input$region_range)
    len   <- .pt_chr_length(chr_lengths, input$region_chr, Inf)
    start <- as.integer(input$region_range[1])
    end   <- as.integer(input$region_range[2])
    if (is.na(start) || start < 1)
      return(list(valid = FALSE, msg = "Start must be \u2265 1"))
    if (is.na(end) || end <= start)
      return(list(valid = FALSE, msg = "End must be > start"))
    if (end > len)
      return(list(valid = FALSE, msg = "End exceeds chromosome length"))
    list(valid = TRUE, chrom = input$region_chr, start = start, end = end)
  })

  output$region_error_ui <- renderUI({
    if (is.null(input$region_mode) || input$region_mode != "coords")
      return(NULL)
    parse_errs <- c(state$start_parse_err, state$end_parse_err)
    if (length(parse_errs) > 0)
      return(div(class = "alert alert-warning",
                 paste(parse_errs, collapse = " | ")))
    v <- tryCatch(region_valid(), error = function(e) NULL)
    if (!is.null(v) && !isTRUE(v$valid))
      div(class = "alert alert-danger", v$msg)
  })

  list(chr_len_r = chr_len_r, region_valid = region_valid)
}

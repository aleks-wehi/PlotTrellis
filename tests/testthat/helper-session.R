# A mock session that records what the server sends the page: custom messages (in
# `$custom`, as list(type, message)) and input updates (in `$input`, as list(id, message)).
recording_session <- function() {
  rec <- new.env()
  rec$custom <- list()
  rec$input <- list()
  session <- shiny::MockShinySession$new()
  session$sendCustomMessage <- function(type, message)
    rec$custom[[length(rec$custom) + 1]] <- list(type = type, message = message)
  session$sendInputMessage <- function(inputId, message)
    rec$input[[length(rec$input) + 1]] <- list(id = inputId, message = message)
  rec$session <- session
  rec
}

# The HTML of a renderUI() output read inside testServer(), as one string.
ui_html <- function(x) paste(as.character(x$html), collapse = "")

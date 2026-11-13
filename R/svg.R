# Rewrites an svglite file in place so it edits cleanly in a vector
# editor (Inkscape, Illustrator) and renders the same everywhere. Run on every SVG
# PlotTrellis writes (the app, pt_save(), the downloaded R script), after the device
# is closed. svglite writes one element per line, which the line-wise edits rely on;
# only each line's opening tag is edited, never a label's text.
#
# unlock_text = FALSE keeps svglite's textLength width lock: the app's own display
# uses it so text keeps R's widths in a browser without the figure's font (multi-part
# headings are placed segment by segment). Files people take away unlock it, so an
# edited label grows instead of being squeezed to the old width.
.pt_tidy_svg <- function(path, unlock_text = TRUE) {
  if (!file.exists(path)) return(invisible(path))
  svg <- readLines(path, warn = FALSE)
  # svglite's white page background, and shapes that draw nothing (no fill, stroke:
  # none): plotgardener's page and per-track viewport boxes, and empty Hi-C bins. An
  # editor's "fit page to drawing" snaps to them, so whitespace can't be trimmed.
  # Clip-path rects have no style and stay.
  bg <- grepl("^<rect width='100%' height='100%'", svg)
  shape <- grepl("^<(rect|polygon|polyline|path|circle|line) [^>]*style='", svg)
  style <- ifelse(shape, sub("^[^>]*style='([^']*)'.*$", "\\1", svg), NA)
  draws_nothing <- !is.na(style) & grepl("stroke: none", style, fixed = TRUE) &
    (!grepl("fill:", style, fixed = TRUE) | grepl("fill: none", style, fixed = TRUE))
  svg <- svg[!(bg | draws_nothing)]
  tag <- regexpr("^<[^>]*>", svg)
  open_tags <- regmatches(svg, tag)
  if (unlock_text) open_tags <- gsub(" textLength='[^']*' lengthAdjust='[^']*'", "", open_tags)
  # A generic fallback for viewers without the font; every font choice is a sans.
  open_tags <- gsub("font-family: ([^;']*[^;' ]);", "font-family: \\1, sans-serif;", open_tags)
  open_tags <- gsub(", sans-serif, sans-serif;", ", sans-serif;", open_tags, fixed = TRUE)
  regmatches(svg, tag) <- open_tags
  # Keeps the leading space of a heading's later segments in renderers that ignore
  # svglite's CSS white-space rule (rsvg drops it without this).
  if (!any(grepl("xml:space=", svg, fixed = TRUE))) svg <- sub("^<svg ", "<svg xml:space='preserve' ", svg)
  writeLines(svg, path)
  invisible(path)
}

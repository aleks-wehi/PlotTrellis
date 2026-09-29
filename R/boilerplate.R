`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

SCALE_BAR_TYPES <- c("bw", "hic", "bw_diff_stacked")

# Flat constant, not fontsize-derived -- matches the pre-refactor inline axis math.
GENOME_AXIS_RESERVE <- 0.35

nice_bp <- function(r) {
  m    <- 10^floor(log10(r * 0.15))
  cand <- m * c(1, 2, 5, 10)
  cand[which.min(abs(cand / r - 0.15))]
}
.fmt_num <- function(x, digits = 2) {
  s <- formatC(x, format = "f", digits = digits)
  s <- sub("0+$", "", s); sub("\\.$", "", s)
}
format_bp_generic <- function(n, digits = 2, space_before_bp = FALSE) {
  if (n >= 1e6) paste0(.fmt_num(n / 1e6, digits), "Mb")
  else if (n >= 1e3) paste0(.fmt_num(n / 1e3, digits), "Kb")
  else paste0(n, if (space_before_bp) " bp" else "bp")
}

.gene_missing <- function(gene) is.null(gene) || is.na(gene)

# transcript_filter.R's functions expect ORGDB_PKG/TX_NAME_LOOKUP_PATH/txdb/
# MIN_DISPLAY_TPM as bare globals, not parameters -- call this before any of them.
.set_transcript_filter_globals <- function(settings) {
  assign("ORGDB_PKG",           settings$orgdb_pkg,                  envir = .GlobalEnv)
  assign("TX_NAME_LOOKUP_PATH", settings$tx_name_lookup_path,        envir = .GlobalEnv)
  assign("txdb",                AnnotationDbi::loadDb(settings$txdb_path), envir = .GlobalEnv)
  assign("MIN_DISPLAY_TPM",     settings$min_display_tpm %||% 1,     envir = .GlobalEnv)
}

# wizard.R's wizard_ui()/wizard_server() expect SETTINGS_JSON/CFG/
# DEFAULT_SETTINGS as bare globals, not parameters (bug 47) -- call this
# before shinyApp(ui = wizard_ui(), server = wizard_server) in app.R.
.set_wizard_globals <- function(settings_json, cfg, default_settings) {
  assign("SETTINGS_JSON",    settings_json,    envir = .GlobalEnv)
  assign("CFG",              cfg,              envir = .GlobalEnv)
  assign("DEFAULT_SETTINGS", default_settings, envir = .GlobalEnv)
}

mixed_text <- function(segments, x, y, fontsize, fontfamily) {
  seg_w <- vapply(segments, function(seg) {
    fface <- if (is.null(seg$fontface)) "plain" else seg$fontface
    grid::convertWidth(grid::grobWidth(grid::textGrob(seg$text, gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily, fontface = fface))), "inches", valueOnly = TRUE)
  }, numeric(1))
  text_h <- grid::convertHeight(grid::grobHeight(grid::textGrob("Xg", gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily))), "inches", valueOnly = TRUE)
  pad <- 0.03
  plotRect(x = x - pad, y = y - pad, width = sum(seg_w) + 2 * pad, height = text_h + 2 * pad,
           just = c("left", "top"), default.units = "inches", fill = "white", linecolor = NA)
  cur_x <- x
  for (si in seq_along(segments)) {
    seg   <- segments[[si]]
    fface <- if (is.null(seg$fontface)) "plain" else seg$fontface
    plotText(label = seg$text, x = cur_x, y = y, just = c("left", "top"),
             fontsize = fontsize, fontfamily = fontfamily, fontface = fface,
             fontcolor = seg$colour)
    cur_x <- cur_x + seg_w[si]
  }
}

draw_value_axis <- function(range, x, y_top, h, fontsize, fontfamily, fontcolor, integer_ticks = FALSE) {
  if (diff(range) <= 0) return(invisible(NULL))
  ticks <- pretty(range)
  if (integer_ticks) ticks <- unique(round(ticks))
  ticks <- ticks[ticks >= range[1] & ticks <= range[2]]
  if (length(ticks) == 0) return(invisible(NULL))
  ys <- y_top + h * (1 - (ticks - range[1]) / diff(range))
  plotSegments(x, y_top, x, y_top + h, default.units = "inches", lwd = 1, linecolor = fontcolor)
  for (i in seq_along(ticks)) {
    plotSegments(x - 0.04, ys[i], x, ys[i], default.units = "inches", lwd = 1, linecolor = fontcolor)
    plotText(label = ticks[i], x = x - 0.06, y = ys[i], just = c("right", "center"), fontsize = fontsize, fontfamily = fontfamily, fontcolor = fontcolor)
  }
}

no_data_text <- function(x, y, w, h, fontsize, fontfamily, fontcolor, label = "No data") {
  plotText(label = label, x = x + w / 2, y = y + h / 2, just = "center",
           fontsize = fontsize, fontfamily = fontfamily, fontcolor = fontcolor)
}

# Every entry is a GROUP of equivalent spellings for one chromosome; a spelling
# with a capture group (e.g. "chr(\\d+)" paired with "(\\d+)") means the two
# spellings agree on the same digit, so resolve_chrom() only needs one group
# entry to cover every digit/X/Y chromosome instead of enumerating chr1..chr22.
CHROM_ALIASES <- list(
  c("chr([0-9XY]+)", "([0-9XY]+)"),
  c("chrM", "MT", "M", "chrMT")
)

# Resolves region_params$chrom against `vocab` (a file's own chromosome
# names). Returns vocab's own spelling, or NULL (soft "no data", never a
# crash) if nothing matches. vocab = NULL trusts the caller's chrom as-is.
resolve_chrom <- function(chrom, vocab) {
  if (is.null(vocab)) return(chrom)
  match_in_vocab <- function(cand) {
    hit <- vocab[tolower(vocab) == tolower(cand)]
    if (length(hit) > 0) hit[1] else NULL
  }
  direct <- match_in_vocab(chrom)
  if (!is.null(direct)) return(direct)

  for (group in CHROM_ALIASES) {
    for (spelling in group) {
      m <- regmatches(chrom, regexec(paste0("^", spelling, "$"), chrom, ignore.case = TRUE))[[1]]
      if (length(m) == 0) next
      captured <- if (length(m) > 1) m[2] else NA_character_
      for (other in setdiff(group, spelling)) {
        has_capture <- grepl("(", other, fixed = TRUE)
        cand <- if (!has_capture) other else if (!is.na(captured)) sub("\\([^)]*\\)", captured, other) else NA_character_
        if (!is.na(cand)) {
          hit <- match_in_vocab(cand)
          if (!is.null(hit)) return(hit)
        }
      }
    }
  }
  NULL
}

draw_scale_bar <- function(ctx) {
  sc <- ctx$scale
  by <- ctx$y + 0.12
  plotSegments(sc$bx0, by,        sc$bx1, by,        default.units = "inches", lwd = 1.5, linecolor = ctx$fontcol)
  plotSegments(sc$bx0, by - 0.05, sc$bx0, by + 0.05, default.units = "inches", lwd = 1.5, linecolor = ctx$fontcol)
  plotSegments(sc$bx1, by - 0.05, sc$bx1, by + 0.05, default.units = "inches", lwd = 1.5, linecolor = ctx$fontcol)
  plotText(sc$label, x = sc$x_mid, y = by - 0.02, just = c("center", "bottom"),
           fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol)
}

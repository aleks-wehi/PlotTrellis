settings_protein_domains <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- 2.0
  list(dict = dict, entry = entry)
}

# No load_protein_domains() -- height is fixed regardless of the UniProt
# lookup outcome, so the whole lookup happens here rather than in settings.
render_protein_domains <- function(entry, dict, ctx) {
  gene <- ctx$region_params$gene
  res <- if (.gene_missing(gene)) {
    list(status = "no_gene")
  } else {
    gene_protein_domains(gene, ctx$txdb, ctx$gtf_path, ctx$orgdb_pkg, ctx$genome_label)
  }

  plotText(label = "Canonical protein annotation", x = ctx$track_x + 0.05, y = ctx$y + 0.02,
           just = c("left", "top"), fontsize = ctx$fontsize - 1, fontface = "bold",
           fontfamily = ctx$font, fontcolor = ctx$fontcol)
  header_h <- 0.2

  pd_text_w <- function(label, fontsize, fontfamily) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily))), "inches", valueOnly = TRUE)
  }
  pd_truncate_label <- function(note, suffix, maxw, fontsize, fontfamily) {
    full <- paste0(note, suffix)
    if (pd_text_w(full, fontsize, fontfamily) <= maxw) return(full)
    n <- nchar(note)
    repeat {
      candidate <- paste0(substr(note, 1, n), "\u2026", suffix)
      if (n <= 1 || pd_text_w(candidate, fontsize, fontfamily) <= maxw) return(candidate)
      n <- n - 1
    }
  }

  if (res$status != "ok") {
    msg <- switch(res$status,
      no_gene = "No gene selected for this region",
      no_uniprot_entry = "No reviewed UniProt entry found for this gene",
      uniprot_unreachable = "UniProt could not be reached or returned an error (no internet on this machine?) -- retries automatically after 2 minutes",
      no_domains = "No annotated domain/region/motif features for this gene's canonical isoform",
      no_canonical_transcript = "UniProt has domain/region/motif features for this gene, but neither the CCDS ID nor a protein-length match could resolve them to a transcript",
      "Protein domain lookup failed (network issue?) -- will retry on next render")
    plotText(label = msg, x = ctx$track_x + ctx$track_w / 2, y = ctx$y + header_h + (ctx$height - header_h) / 2,
             just = "center", fontsize = ctx$fontsize - 1, fontfamily = ctx$font, fontcolor = ctx$fontcol)
    return(invisible(NULL))
  }

  method_label <- c(
    ccds_id_join = "via CCDS ID",
    length_match_fallback = "via protein-length match (CCDS not tagged in GTF yet)",
    length_match_fallback_ambiguous = "via protein-length match, ambiguous -- multiple transcripts share this length, first taken"
  )
  subtitle <- sprintf("canonical: %s (%s)", res$canonical_tx, method_label[[res$method]])
  plotText(label = subtitle, x = ctx$track_x + 0.05, y = ctx$y + 0.15, just = c("left", "top"),
           fontsize = ctx$fontsize - 3, fontfamily = ctx$font, fontcolor = "grey45")

  domains <- res$domains
  dids <- unique(domains$domain_id)
  dids <- dids[order(vapply(dids, function(d) min(domains$start[domains$domain_id == d]), numeric(1)))]
  n <- length(dids)

  if (n == 0) {
    plotText(label = "No annotated domain/region/motif features for this gene's canonical isoform",
             x = ctx$track_x + ctx$track_w / 2, y = ctx$y + (ctx$height + header_h) / 2, just = "center",
             fontsize = ctx$fontsize - 1, fontfamily = ctx$font, fontcolor = ctx$fontcol)
    return(invisible(NULL))
  }

  header_h <- 0.33  # grows to fit the type legend below the subtitle
  present_types <- intersect(names(PD_TYPE_COLOURS), unique(domains$type))
  legend_x <- ctx$track_x + 0.05
  for (ty in present_types) {
    plotRect(x = legend_x, y = ctx$y + 0.27, width = 0.09, height = 0.09, just = c("left", "center"),
             default.units = "inches", fill = PD_TYPE_COLOURS[[ty]], linecolor = NA)
    legend_x <- legend_x + 0.09 + 0.04
    plotText(label = ty, x = legend_x, y = ctx$y + 0.27, just = c("left", "center"),
             fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)
    legend_x <- legend_x + pd_text_w(ty, ctx$fontsize - 2, ctx$font) + 0.14
  }

  # PAGE_W - TRACK_X - TRACK_W is always exactly 1.0in (render_all's fixed right margin).
  label_maxw <- 1.0 - 0.05 - 0.05
  row_h <- (ctx$height - header_h) / n

  for (k in seq_len(n)) {
    did  <- dids[k]
    seg  <- domains[domains$domain_id == did, ]
    seg  <- seg[order(seg$start), ]
    ty   <- seg$type[1]
    col  <- PD_TYPE_COLOURS[[ty]]
    nseg <- nrow(seg)
    y0   <- ctx$y + header_h + (k - 1) * row_h
    ymid <- y0 + row_h / 2

    if (nseg > 1) {
      x0 <- ctx$track_x + (min(seg$start) - ctx$region_params$start) / (ctx$region_params$end - ctx$region_params$start) * ctx$track_w
      x1 <- ctx$track_x + (max(seg$end)   - ctx$region_params$start) / (ctx$region_params$end - ctx$region_params$start) * ctx$track_w
      plotSegments(x0 = x0, x1 = x1, y0 = ymid, y1 = ymid, default.units = "inches", linecolor = col, lwd = 1, lty = 2)
      plotText(label = sprintf("%d exons", nseg), x = ctx$track_x - 0.05, y = ymid, just = c("right", "center"),
               fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = "grey50")
    }

    for (s in seq_len(nrow(seg))) {
      sx0 <- ctx$track_x + (seg$start[s] - ctx$region_params$start) / (ctx$region_params$end - ctx$region_params$start) * ctx$track_w
      sx1 <- ctx$track_x + (seg$end[s]   - ctx$region_params$start) / (ctx$region_params$end - ctx$region_params$start) * ctx$track_w
      w  <- max(sx1 - sx0, 0.012)
      xc <- (sx0 + sx1) / 2
      plotRect(x = xc - w / 2, y = ymid - 0.055, width = w, height = 0.11, just = c("left", "top"),
               default.units = "inches", fill = col, linecolor = NA)
    }

    label <- pd_truncate_label(seg$note[1], "", label_maxw, ctx$fontsize - 1.5, ctx$font)
    plotText(label = label, x = ctx$track_x + ctx$track_w + 0.05, y = ymid, just = c("left", "center"),
             fontsize = ctx$fontsize - 1.5, fontfamily = ctx$font, fontcolor = col)
  }
  invisible(NULL)
}

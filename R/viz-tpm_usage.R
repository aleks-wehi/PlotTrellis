# dict["tpm_range"] is a shared, figure-wide axis max (not per-instance); unlike production, this pipeline only ever sees one entry at a time, so it doesn't merge contiguous tpm_usage tracks into one panel -- correct for a lone track (the only case tested), a real gap for 2+ adjacent ones (BUGS.md #46).
settings_tpm_usage <- function(dict, entry, settings, track_w, region_params) {
  entry <- load_transcripts(entry, region_params, settings)
  gene <- region_params$gene
  if (.gene_missing(gene)) {
    entry$data$all <- data.frame(transcript_id = character(0), transcript_name = character(0),
                                  max_tpm = numeric(0), mean_tpm = numeric(0), proportion = numeric(0))
    entry$data$total_tpm <- 0
  } else {
    tpm_samples <- list(list(sample_id = entry$sample, junction_path = entry$junction_path, quant_path = entry$quant_path))
    all_tx <- filtered_transcripts(gene = gene, samples = tpm_samples, gtf_path = settings$gtf_path)
    entry$data$all <- all_tx
    entry$data$total_tpm <- attr(all_tx, "total_tpm")
  }

  if (is.null(entry$height)) entry$height <- 1.35  # run_len<=1 branch -- see scope-limitation note above
  cur <- dict[["tpm_range"]]
  this_total <- entry$data$total_tpm %||% 0
  dict[["tpm_range"]] <- if (is.null(cur)) this_total else max(cur, this_total)
  list(dict = dict, entry = entry)
}

render_tpm_usage <- function(entry, dict, ctx) {
  TPM_AXIS_H <- 0.25

  gene <- ctx$region_params$gene
  if (.gene_missing(gene)) {
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol)
    return(invisible(NULL))
  }

  named_union <- entry$data$named
  n_excluded  <- entry$data$n_annotated - length(named_union)
  has_other   <- n_excluded > 0
  other_label <- if (has_other) sprintf("Other (n=%d)", n_excluded) else NA_character_
  colours <- transcript_colours(named_union)
  if (has_other) colours <- c(colours, setNames(TX_OTHER_COLOUR, other_label))

  all_tx <- entry$data$all
  tpm_named  <- all_tx[all_tx$transcript_name %in% named_union, ]
  other_tpm  <- max(0, entry$data$total_tpm - sum(tpm_named$mean_tpm))
  other_prop <- if (entry$data$total_tpm > 0) other_tpm / entry$data$total_tpm else 0

  n_rows_needed <- nrow(tpm_named) + as.integer(has_other)
  combined <- data.frame(
    sample_name = rep(entry$sample, n_rows_needed),
    transcript_name = factor(if (has_other) c(tpm_named$transcript_name, other_label) else tpm_named$transcript_name,
                              levels = if (has_other) c(named_union, other_label) else named_union),
    tpm = if (has_other) c(tpm_named$mean_tpm, other_tpm) else tpm_named$mean_tpm,
    proportion = if (has_other) c(tpm_named$proportion, other_prop) else tpm_named$proportion
  )
  combined$sample_name <- factor(combined$sample_name, levels = c(entry$sample))

  tpm_range <- max(1, dict[["tpm_range"]]) * (10 / 9)

  bar_h <- max(0.15, ctx$height - TPM_AXIS_H)
  row_center <- ctx$y + bar_h / 2

  combined$tpm_label <- ifelse(combined$proportion >= 0.05, sprintf("%.0f%%", combined$proportion * 100), "")

  # Text sizes scale with the fontsize setting; equal to the old fixed sizes at the default 8.
  gg <- ggplot2::ggplot(combined, ggplot2::aes(x = sample_name, y = tpm, fill = transcript_name)) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.6, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = tpm_label),
      position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE), colour = "white", size = 2.6 * ctx$fontsize / 8, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colours, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, tpm_range), expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::theme_minimal(base_size = ctx$fontsize) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 0, 0),
      axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right", legend.key.size = grid::unit(0.25, "cm"), legend.text = ggplot2::element_text(size = ctx$fontsize - 1))

  label_w <- 1.4
  plotGG(plot = gg, x = ctx$track_x + label_w, y = ctx$y, width = ctx$track_w - label_w, height = bar_h, just = c("left", "top"))

  plotText(label = "Transcript counts", x = ctx$track_x + 0.05, y = ctx$y + 0.02, just = c("left", "top"),
           fontsize = ctx$fontsize - 1, fontface = "bold", fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText(label = entry$sample, x = ctx$track_x + 0.05, y = row_center, just = c("left", "center"),
           fontsize = ctx$fontsize - 1, fontfamily = ctx$font, fontcolor = entry$colour)

  axis_y <- ctx$y + bar_h
  ticks <- pretty(c(0, tpm_range))
  ticks <- ticks[ticks >= 0 & ticks <= tpm_range]
  tick_x <- (ctx$track_x + label_w) + (ticks / tpm_range) * (ctx$track_w - label_w)
  for (k in seq_along(ticks)) {
    plotSegments(x0 = tick_x[k], y0 = axis_y, x1 = tick_x[k], y1 = axis_y + 0.04, default.units = "inches", linecolor = "grey50", lwd = 0.5)
    plotText(label = ticks[k], x = tick_x[k], y = axis_y + 0.05, just = c("center", "top"),
             fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  }
  plotText(label = "TPM", x = ctx$track_x + label_w + (ctx$track_w - label_w) / 2, y = ctx$y + ctx$height,
           just = c("center", "bottom"), fontsize = ctx$fontsize - 1.5, fontfamily = ctx$font, fontcolor = ctx$fontcol)

  invisible(NULL)
}

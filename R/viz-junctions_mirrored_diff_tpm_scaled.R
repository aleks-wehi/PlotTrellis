# Fixed height, no dict fold -- the pool DOES happen here (not in render_<viz>)
# because height_scale_treat needs filtered_transcripts(), which needs
# .set_transcript_filter_globals(settings); render_<viz> never receives settings.
settings_junctions_mirrored_diff_tpm_scaled <- function(dict, entry, settings, track_w, region_params) {
  gene    <- region_params$gene
  no_gene <- is.null(gene) || is.na(gene) || nchar(trimws(gene)) == 0

  raw_ctrl  <- load_junctions(entry$control_path,   region_params)
  raw_treat <- load_junctions(entry$treatment_path, region_params)

  if (no_gene) {
    entry$data <- list(no_gene = TRUE)
  } else {
    .set_transcript_filter_globals(settings)

    ctrl_tpm  <- attr(filtered_transcripts(gene = gene, samples = list(list(sample_id = entry$control_label,   quant_path = entry$control_quant)),   gtf_path = settings$gtf_path), "total_tpm")
    treat_tpm <- attr(filtered_transcripts(gene = gene, samples = list(list(sample_id = entry$treatment_label, quant_path = entry$treatment_quant)), gtf_path = settings$gtf_path), "total_tpm")
    # Reference is the anchor (scale 1); falls back to unscaled if Comparison's gene TPM is 0.
    height_scale_treat <- if (treat_tpm > 0) ctrl_tpm / treat_tpm else 1

    entry$data <- list(no_gene = FALSE, raw_ctrl = raw_ctrl, raw_treat = raw_treat, height_scale_treat = height_scale_treat)
  }

  if (is.null(entry$height)) entry$height <- 1.15
  list(dict = dict, entry = entry)
}

render_junctions_mirrored_diff_tpm_scaled <- function(entry, dict, ctx) {
  d <- entry$data
  if (isTRUE(d$no_gene)) {
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol)
    return(invisible(NULL))
  }

  JUNCMDT_PSEUDOCOUNT <- 1
  JUNCMDT_LEGEND_N   <- 100
  JUNCMDT_LEGEND_IMG <- matrix(rev(sapply(seq(0, 1, length.out = JUNCMDT_LEGEND_N), function(a) grDevices::adjustcolor("black", alpha.f = a))), ncol = 1)

  # Differs from junctions_mirrored_diff's juncmd_to_arches() in one respect:
  # `height_scale` multiplies nSupport (1 for Reference, ctrl/treat TPM ratio
  # for Comparison) before it becomes archHeight; opacity is looked up, not recomputed.
  juncmdt_to_arches <- function(chrom, pooled, merged, cap, height_scale) {
    if (nrow(pooled) == 0) {
      return(data.frame(chrom1 = character(0), start1 = integer(0), end1 = integer(0),
                         chrom2 = character(0), start2 = integer(0), end2 = integer(0),
                         nSupport = integer(0), opacity = numeric(0)))
    }
    idx <- match(paste(pooled$intron_start, pooled$intron_end),
                 paste(merged$intron_start, merged$intron_end))
    opacity <- pmin(merged$abs_lfc[idx] / cap, 1)
    data.frame(
      chrom1 = chrom, start1 = pooled$intron_start, end1 = pooled$intron_start + 1L,
      chrom2 = chrom, start2 = pooled$intron_end,   end2 = pooled$intron_end + 1L,
      nSupport = pooled$nSupport * height_scale,
      opacity  = opacity
    )
  }

  raw_ctrl  <- d$raw_ctrl
  raw_treat <- d$raw_treat
  min_support  <- ctx$min_junction_support
  pooled_ctrl  <- raw_ctrl[raw_ctrl$nSupport   >= min_support, ]
  pooled_treat <- raw_treat[raw_treat$nSupport >= min_support, ]
  cap <- ctx$junction_opacity_cap

  merged <- merge(raw_treat, raw_ctrl, by = c("intron_start", "intron_end"), all = TRUE, suffixes = c("_treat", "_ctrl"))
  merged$nSupport_treat[is.na(merged$nSupport_treat)] <- 0
  merged$nSupport_ctrl[is.na(merged$nSupport_ctrl)]   <- 0
  # Opacity uses Comparison's TPM-SCALED count, so it agrees with archHeight on
  # how "TPM-corrected" this pair's Comparison count is.
  merged$nSupport_treat_scaled <- merged$nSupport_treat * d$height_scale_treat
  merged$abs_lfc <- abs(log2((merged$nSupport_treat_scaled + JUNCMDT_PSEUDOCOUNT) / (merged$nSupport_ctrl + JUNCMDT_PSEUDOCOUNT)))

  arches_treat <- juncmdt_to_arches(ctx$region_params$chrom, pooled_treat, merged, cap, d$height_scale_treat)
  arches_ctrl  <- juncmdt_to_arches(ctx$region_params$chrom, pooled_ctrl,  merged, cap, 1)
  shared_max <- if (nrow(arches_treat) > 0 || nrow(arches_ctrl) > 0) max(arches_treat$nSupport, arches_ctrl$nSupport) else 1

  d <- list(arches_treat = arches_treat, arches_ctrl = arches_ctrl, shared_max = shared_max, cap = cap,
            height_scale_treat = d$height_scale_treat)

  header_h    <- 0.3
  content_top <- ctx$y + header_h
  half_h      <- (ctx$height - header_h) / 2
  mid_y       <- content_top + half_h

  treat_fill <- vapply(d$arches_treat$opacity, function(a) grDevices::adjustcolor(entry$treatment_colour, alpha.f = a), character(1))
  ctrl_fill  <- vapply(d$arches_ctrl$opacity,  function(a) grDevices::adjustcolor(entry$control_colour,   alpha.f = a), character(1))

  if (nrow(d$arches_treat) > 0) {
    plotPairsArches(
      data = d$arches_treat, params = ctx$region,
      archHeight = "nSupport", range = c(0, d$shared_max), flip = FALSE,
      fill = treat_fill, linecolor = treat_fill, alpha = 1, baseline = TRUE,
      x = ctx$track_x, y = content_top, width = ctx$track_w, height = half_h
    )
  }
  if (nrow(d$arches_ctrl) > 0) {
    plotPairsArches(
      data = d$arches_ctrl, params = ctx$region,
      archHeight = "nSupport", range = c(0, d$shared_max), flip = TRUE,
      fill = ctrl_fill, linecolor = ctrl_fill, alpha = 1, baseline = TRUE,
      x = ctx$track_x, y = mid_y, width = ctx$track_w, height = half_h
    )
  }

  mixed_text(
    list(
      list(text = paste0(entry$name, " mirrored diff (TPM-scaled)"), colour = ctx$fontcol, fontface = "bold"),
      list(text = " · ", colour = ctx$fontcol),
      list(text = entry$treatment_label, colour = entry$treatment_colour),
      list(text = sprintf(" (Cmp, ×%.2f gene-TPM scaled) vs ", d$height_scale_treat), colour = ctx$fontcol),
      list(text = entry$control_label, colour = entry$control_colour),
      list(text = " (Ref)", colour = ctx$fontcol)
    ),
    x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
  )

  axis_x         <- ctx$track_x - 0.08
  content_bottom <- mid_y + half_h
  plotSegments(axis_x, content_top, axis_x, content_bottom, default.units = "inches", lwd = 1, linecolor = ctx$fontcol)
  plotSegments(axis_x - 0.04, content_top,    axis_x, content_top,    default.units = "inches", lwd = 1, linecolor = ctx$fontcol)
  plotSegments(axis_x - 0.04, mid_y,          axis_x, mid_y,          default.units = "inches", lwd = 1, linecolor = ctx$fontcol)
  plotSegments(axis_x - 0.04, content_bottom, axis_x, content_bottom, default.units = "inches", lwd = 1, linecolor = ctx$fontcol)
  plotText(sprintf("%d", as.integer(round(d$shared_max))), x = axis_x - 0.06, y = content_top,    just = c("right", "center"), fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText("0",                                             x = axis_x - 0.06, y = mid_y,          just = c("right", "center"), fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText(sprintf("%d", as.integer(round(d$shared_max))), x = axis_x - 0.06, y = content_bottom, just = c("right", "center"), fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)

  legend_h <- min(ctx$height * 0.6, 1.0)
  legend_x <- ctx$track_x + ctx$track_w + 0.2
  legend_y <- ctx$y + (ctx$height - legend_h) / 2
  plotText("Opacity", x = legend_x + 0.09, y = legend_y - 0.08, just = c("center", "bottom"), fontsize = ctx$fontsize - 1, fontfamily = ctx$font, fontcolor = ctx$fontcol, fontface = "bold")
  plotRaster(image = JUNCMDT_LEGEND_IMG, x = legend_x, y = legend_y, width = 0.18, height = legend_h, just = c("left", "top"), default.units = "inches", interpolate = TRUE)
  plotText(sprintf("%.2f log2FC", d$cap), x = legend_x + 0.22, y = legend_y, just = c("left", "top"), fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText("(full)", x = legend_x + 0.22, y = legend_y + 0.11, just = c("left", "top"), fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText("0", x = legend_x + 0.22, y = legend_y + legend_h, just = c("left", "bottom"), fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = ctx$fontcol)

  invisible(NULL)
}

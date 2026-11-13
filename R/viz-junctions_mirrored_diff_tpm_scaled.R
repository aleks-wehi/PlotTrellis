info_junctions_mirrored_diff_tpm_scaled <- function() list(
  kind = "comparison", data_type = "paired-end", order = 30, label = "{name} mirrored diff (TPM-scaled)",
  api = "{name}_diff_tpm", filter = c(comparison = "Comparisons"), needs_gene = TRUE
)

# Offered only when both samples have a quantification, whose paths its items carry.
enrich_junctions_mirrored_diff_tpm_scaled <- function(items, control_row, treatment_row) {
  ctrl_quant  <- sample_data_path(control_row, "quantdir")
  treat_quant <- sample_data_path(treatment_row, "quantdir")
  if (is.na(ctrl_quant) || is.na(treat_quant)) return(list())
  lapply(items, function(it) modifyList(it, list(control_quant = ctrl_quant, treatment_quant = treat_quant)))
}

# The mirrored junction comparison with Comparison's read counts scaled by the two samples'
# gene-TPM ratio. The junctions are loaded and the scale computed here, not in render_<viz>,
# because height_scale_treat needs filtered_transcripts(), which needs the annotation
# from .pt_annotation(settings); render_<viz> never receives settings.
settings_junctions_mirrored_diff_tpm_scaled <- function(dict, entry, settings, track_w, region_params) {
  gene    <- region_params$gene
  no_gene <- .pt_no_gene(gene)

  raw_ctrl  <- load_junctions(entry$control_path,   region_params)
  raw_treat <- load_junctions(entry$treatment_path, region_params)

  if (no_gene) {
    entry$data <- list(no_gene = TRUE)
  } else {
    anno <- .pt_annotation(settings)

    ctrl  <- filtered_transcripts(gene = gene, samples = list(list(sample_id = entry$control_label,   quant_path = entry$control_quant)),   anno = anno)
    treat <- filtered_transcripts(gene = gene, samples = list(list(sample_id = entry$treatment_label, quant_path = entry$treatment_quant)), anno = anno)
    bad <- c(attr(ctrl, "tpm_not_numeric"), attr(treat, "tpm_not_numeric"))
    if (length(bad)) stop("the TPM column of a quant file is not all numbers (corrupt file?): ", bad[1], call. = FALSE)
    ctrl_tpm  <- attr(ctrl, "total_tpm")
    treat_tpm <- attr(treat, "total_tpm")
    # Reference is the anchor (scale 1); falls back to unscaled if Comparison's gene TPM is 0.
    height_scale_treat <- if (treat_tpm > 0) ctrl_tpm / treat_tpm else 1

    entry$data <- list(no_gene = FALSE, raw_ctrl = raw_ctrl, raw_treat = raw_treat, height_scale_treat = height_scale_treat)
  }

  if (is.null(entry$height)) entry$height <- JUNCTION_HEADING_H + 2 * settings$junction_height
  list(dict = dict, entry = entry)
}

render_junctions_mirrored_diff_tpm_scaled <- function(entry, dict, ctx) {
  d <- entry$data
  if (isTRUE(d$no_gene)) {
    .pt_status_text(ctx, .pd_status_text("no_gene"))
    return(invisible(NULL))
  }
  # Comparison's read counts times the Reference/Comparison gene-TPM ratio, in the arc heights
  # and in the fold change behind the opacity alike.
  arches <- mirrored_junction_arches(
    treat = d$raw_treat, ctrl = d$raw_ctrl, chrom = ctx$region_params$chrom,
    min_support = ctx$min_junction_support, cap = ctx$junction_opacity_cap, treat_scale = d$height_scale_treat
  )
  render_mirrored_junction_track(arches, " mirrored diff (TPM-scaled)",
                                 sprintf(" (Cmp, \u00D7%.2f gene-TPM scaled) vs ", d$height_scale_treat), entry, ctx)
}

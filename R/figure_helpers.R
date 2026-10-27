# Figure-building helpers shared by the composer app (R/app.R) and the code API
# (R/api.R), so a figure made either way gets the same entries, settings and
# colours. They used to live inside .pt_composer_app()'s closure, bound to its
# GENOME/TRACK_ITEMS/ORGDB_PKG globals; now those are arguments.

PALETTE <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a",
             "#a6761d","#e6ab02","#66a61e","#e7298a","#7570b3")

# Treat NULL/NA/non-finite/list alike as "unset" (list: jsonlite::toJSON() used to round-trip a NULL setting as {}, fixed at the write site in wizard.R's wiz_save, but a settings.json written before that fix still reads {} back as an empty list via fromJSON()); shared by RENDER_CFG's zrange_cap and hic_height, both falling back to auto.
clean_optional_numeric <- function(v) {
  if (is.null(v) || is.list(v) || is.na(v) || !is.finite(v)) NULL else v
}

# Every render/lookup call site's single "cfg" object (Phase C.3 decision 7,
# replaces the old cur_settings() reactive + settings_modal()'s session-only
# input$s_... overrides). Bundles a few genome-derived fields
# (orgdb_pkg/genome_label/tx_*_lookup_path) alongside the flat CFG keys purely
# so callers keep getting one object with both.
.pt_render_cfg <- function(CFG, GENOME) list(
  txdb_path   = GENOME$txdb_path,
  gtf_path    = CFG$gtf_path,
  samples     = CFG$samples,
  orgdb_pkg           = GENOME$orgdb,
  genome_label        = GENOME$genome_label,
  tx_name_lookup_path = GENOME$tx_name_lookup_path,
  tx_meta_lookup_path = GENOME$tx_meta_lookup_path,
  data_dir    = CFG$data_dir,
  font        = trimws(CFG$font %||% "DejaVu Sans"),
  fontsize    = CFG$fontsize %||% 8L,
  page_w      = CFG$page_w %||% 9,
  chr_fmt     = CFG$chr_fmt %||% "dot",
  show_track_axes = CFG$show_track_axes %||% TRUE,
  show_title      = CFG$show_title %||% TRUE,
  show_scale_bar  = CFG$show_scale_bar %||% TRUE,
  hic_norm    = trimws(CFG$hic_norm %||% "KR"),
  hic_style   = CFG$hic_style %||% "rectangle",
  zrange_cap  = clean_optional_numeric(CFG$zrange_cap),
  atac_height  = CFG$atac_height %||% 0.85,
  hic_height   = clean_optional_numeric(CFG$hic_height),
  hic_res_mode = CFG$hic_res_mode %||% "auto",
  hic_max_bins_per_view = CFG$hic_max_bins_per_view %||% 500,
  diff_colour_lo  = trimws(CFG$diff_colour_lo %||% "#0033CC"),
  diff_colour_hi  = trimws(CFG$diff_colour_hi %||% "#CC0000"),
  junction_opacity_cap = CFG$junction_opacity_cap %||% 2,
  min_junction_support = CFG$min_junction_support %||% 0,
  min_display_tpm = CFG$min_display_tpm %||% 1,
  reference_colour  = trimws(CFG$reference_colour %||% "#4D4D4D"),
  comparison_colour = trimws(CFG$comparison_colour %||% "#56B4E9"),
  gene_flank_pct  = CFG$gene_flank_pct %||% 15,
  uniprot_reviewed_only = CFG$uniprot_reviewed_only %||% FALSE,
  arc_clip_partial = CFG$arc_clip_partial %||% FALSE,
  junction_height  = CFG$junction_height %||% 0.85,
  diff_shared_alpha = CFG$diff_shared_alpha %||% 0.15
)

# render_all()'s settings, from a .pt_render_cfg() object (the app adds title/output_path per render).
.pt_settings_base <- function(cfg) list(
  page_w    = cfg$page_w,
  font      = cfg$font,
  fontsize  = as.integer(cfg$fontsize),
  hic_norm  = cfg$hic_norm,
  zrange    = if (is.null(cfg$zrange_cap) || is.na(cfg$zrange_cap))
                NULL else c(0, cfg$zrange_cap),
  hic_style = cfg$hic_style,
  hic_res_mode = cfg$hic_res_mode,
  hic_height   = cfg$hic_height,
  hic_max_bins_per_view = cfg$hic_max_bins_per_view,
  atac_height = cfg$atac_height,
  chr_fmt   = cfg$chr_fmt,
  show_track_axes = cfg$show_track_axes,
  show_title      = cfg$show_title,
  show_scale_bar  = cfg$show_scale_bar,
  txdb_path = cfg$txdb_path,
  gtf_path  = cfg$gtf_path,
  orgdb_pkg           = cfg$orgdb_pkg,
  genome_label        = cfg$genome_label,
  tx_name_lookup_path = cfg$tx_name_lookup_path,
  tx_meta_lookup_path = cfg$tx_meta_lookup_path,
  data_dir = cfg$data_dir,
  samples  = cfg$samples,
  diff_colour_lo  = cfg$diff_colour_lo,
  diff_colour_hi  = cfg$diff_colour_hi,
  junction_opacity_cap = cfg$junction_opacity_cap,
  min_junction_support = cfg$min_junction_support,
  min_display_tpm = cfg$min_display_tpm,
  reference_colour  = cfg$reference_colour,
  comparison_colour = cfg$comparison_colour,
  arc_clip_partial  = cfg$arc_clip_partial,
  junction_height   = cfg$junction_height,
  diff_shared_alpha = cfg$diff_shared_alpha
)

# Gene query (symbol in any case, or Ensembl gene ID; see resolve_gene_query()) -> coords,
# the canonical symbol every track/panel keys on, and a user-facing error.
.pt_resolve_gene_region <- function(query, txdb, flank_pct, orgdb_pkg, all_genes, genome_label) {
  q <- resolve_gene_query(query, orgdb_pkg, all_genes)
  result <- if (is.na(q$symbol)) NULL else
    gene_to_coords(q$symbol, txdb = txdb, flank_pct = flank_pct, orgdb_pkg = orgdb_pkg, ensembl_id = q$ensembl)
  if (is.null(result))
    list(region = NULL, symbol = NA_character_,
         error = q$error %||% paste0("'", query, "' not found in ", genome_label,
                                     ": enter a gene symbol (e.g. Pax5) or an Ensembl gene ID"))
  else
    list(region = result, symbol = q$symbol, error = NULL,
         note = if (length(q$also)) paste0(q$ensembl, " is also named ", paste(q$also, collapse = ", ")))
}

# Per-sample colour: the samples.txt colour column when present, else PALETTE in sample order.
.pt_sample_colours <- function(samples_df) {
  samples <- unique(vapply(samples_df, function(r) r$display, character(1)))
  setNames(
    vapply(seq_along(samples), function(i) {
      s   <- samples[i]
      row <- Find(function(r) r$display == s, samples_df)
      if (!is.null(row) && !is.na(row$colour %||% NA_character_)) return(row$colour)
      PALETTE[((i - 1L) %% length(PALETTE)) + 1L]
    }, character(1)),
    samples
  )
}

# Converts TRACK_ITEMS entries (picker selections, utils/track_detection.R's shape) into pipeline/render_all()'s entry shape; height/data are always NULL here since render_all()'s settings_<viz>() functions compute height themselves.
.pt_build_entries <- function(figure_labels, track_items, col_map, settings) {
  Filter(Negate(is.null), lapply(figure_labels, function(lbl) {
    sp <- track_items[[lbl]]
    if (is.null(sp)) return(NULL)

    # Comparison items carry `samples` = c(reference_display, comparison_display)
    # instead of a single colour_key. control_colour/treatment_colour are
    # resolved HERE, from two fixed per-project ROLE colours (settings$
    # reference_colour/comparison_colour), never from either condition's own
    # identity/PALETTE colour -- "which one is the baseline" reads the same
    # way in every comparison track regardless of the samples involved.
    if (!is.null(sp$samples) && length(sp$samples) == 2) {
      known <- c("key", "label", "type", "sample", "file", "colour_key", "samples",
                 "control_label", "treatment_label", "control_path", "treatment_path")
      extra <- sp[setdiff(names(sp), known)]  # e.g. junctions_mirrored_diff_tpm_scaled's control_quant/treatment_quant
      return(c(list(
        viz = sp$type, name = sp$name,
        control_label = sp$control_label, treatment_label = sp$treatment_label,
        control_path  = sp$control_path,  treatment_path  = sp$treatment_path,
        control_colour   = settings$reference_colour  %||% "#4D4D4D",
        treatment_colour = settings$comparison_colour %||% "#56B4E9",
        height = NULL, data = NULL
      ), extra))
    }

    col <- if (!is.na(sp$colour_key %||% NA)) col_map[[sp$colour_key]] else NA_character_
    entry <- list(viz = sp$type, height = NULL, data = NULL)
    if (sp$type == "tpm_usage") {
      entry <- c(entry, list(sample = sp$sample, colour = col, quant_path = sp$quant_path, junction_path = sp$junction_path))
    } else if (sp$type == "sashimi") {
      entry <- c(entry, list(sample = sp$sample, colour = col, name = sp$name, path = sp$junction_path))
    } else if (sp$type %in% c("bw", "hic", "bed")) {
      entry <- c(entry, list(sample = sp$sample, colour = col, name = sp$name, path = sp$file))
      if (sp$type == "hic") entry$hic_style <- settings$hic_style
    }
    entry
  }))
}

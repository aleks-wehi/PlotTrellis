# What the app and the code API offer, read off the modules' own info_<key>()
# declarations (.pt_viz_infos(), see the viz-module-contract in viz-render-all.R).

# A label or API-name template from info_<key>() with `{name}` replaced by the
# data's name.
.pt_fill_name <-
  function(template, name) sub("{name}", name, template, fixed = TRUE)

detect_generic <- function(sample_row, info) {
  instances <- sample_data_instances(sample_row, info$data_type)
  if (length(instances) == 0) return(list())
  safe_s <- gsub("[^A-Za-z0-9._-]", "_", sample_row$display)
  lapply(instances, function(inst) {
    safe_n <- gsub("[^A-Za-z0-9._-]", "_", inst$name)
    item <- list(
      key        = paste0(info$key, "_", safe_n, "_", safe_s),
      label      = paste0(.pt_fill_name(info$label, inst$name), " \u00B7 ",
                          sample_row$display),
      type       = info$key,
      sample     = sample_row$display,
      file       = if (info$path_field == "file") inst$path else NA_character_,
      colour_key = sample_row$display,
      name       = inst$name,
      broken     = !file.exists(inst$path)
    )
    if (info$path_field != "file") item[[info$path_field]] <- inst$path
    item
  })
}

detect_pair_generic <- function(control_row, treatment_row, info) {
  ctrl_insts  <- sample_data_instances(control_row, info$data_type)
  treat_insts <- sample_data_instances(treatment_row, info$data_type)
  ctrl_names  <- vapply(ctrl_insts, `[[`, character(1), "name")
  treat_names <- vapply(treat_insts, `[[`, character(1), "name")
  common <- intersect(ctrl_names, treat_names)
  if (length(common) == 0) return(list())

  safe_c <- gsub("[^A-Za-z0-9._-]", "_", control_row$display)
  safe_t <- gsub("[^A-Za-z0-9._-]", "_", treatment_row$display)

  items <- lapply(common, function(nm) {
    safe_n     <- gsub("[^A-Za-z0-9._-]", "_", nm)
    ctrl_path  <- ctrl_insts[[which(ctrl_names == nm)[1]]]$path
    treat_path <- treat_insts[[which(treat_names == nm)[1]]]$path
    item <- list(
      key = paste0(info$key, "_", safe_n, "_", safe_c, "_vs_", safe_t),
      label = paste0(.pt_fill_name(info$label, nm), " \u00B7 ",
                     treatment_row$display, " (Cmp) vs ",
                     control_row$display, " (Ref)"),
      type = info$key, sample = NA_character_, file = NA_character_,
      colour_key = NA_character_,
      name = nm, samples = c(control_row$display, treatment_row$display),
      control_label = control_row$display,
      treatment_label = treatment_row$display,
      control_path = ctrl_path, treatment_path = treat_path
    )
  })

  enrich <- .pt_module_fn("enrich", info$key)
  if (!is.null(enrich)) items <- enrich(items, control_row, treatment_row)
  items
}

fixed_item <- function(info) {
  list(key = info$key, label = info$label, type = info$key,
       sample = NA_character_, file = NA_character_, colour_key = NA_character_)
}

# Every per-sample track of each sample (sample by sample, types in detection
# order), then the sample-free tracks.
build_track_items <- function(samples_list) {
  infos <- .pt_viz_infos()
  items <- list()
  for (row in samples_list) {
    for (info in Filter(function(info) info$kind == "sample", infos)) {
      detect <- .pt_module_fn("detect", info$key)
      found <- if (is.null(detect))
        detect_generic(row, info) else detect(row, info)
      for (item in found) items[[item$key]] <- item
    }
  }
  for (info in Filter(function(info) info$kind == "fixed", infos)) {
    item <- fixed_item(info)
    items[[item$key]] <- item
  }
  items
}

# Every comparison track of each ordered pair of samples (reference, then
# comparison).
build_comparison_track_items <- function(samples_list) {
  infos <- Filter(function(info) info$kind == "comparison", .pt_viz_infos())
  items <- list()
  n <- length(samples_list)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      control   <- samples_list[[i]]
      treatment <- samples_list[[j]]
      if (identical(control$display, treatment$display)) next
      for (info in infos) {
        for (item in detect_pair_generic(control, treatment, info))
          items[[item$key]] <- item
      }
    }
  }
  items
}

# Track types that draw only for a gene (needs_gene): refused for a region.
.pt_gene_required_types <- function(infos = .pt_viz_infos()) {
  names(Filter(function(info) info$needs_gene, infos))
}

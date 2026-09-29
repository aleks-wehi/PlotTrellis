detect_generic <- function(sample_row, viz, data_type, path_field = "file") {
  instances <- sample_data_instances(sample_row, data_type)
  if (length(instances) == 0) return(list())
  safe_s <- gsub("[^A-Za-z0-9._-]", "_", sample_row$display)
  lapply(instances, function(inst) {
    safe_n <- gsub("[^A-Za-z0-9._-]", "_", inst$name)
    item <- list(
      key        = paste0(viz, "_", safe_n, "_", safe_s),
      label      = paste0(inst$name, " · ", sample_row$display),
      type       = viz,
      sample     = sample_row$display,
      file       = if (path_field == "file") inst$path else NA_character_,
      colour_key = sample_row$display,
      name       = inst$name,
      broken     = !file.exists(inst$path)
    )
    if (path_field != "file") item[[path_field]] <- inst$path
    item
  })
}

detect_pair_generic <- function(control_row, treatment_row, viz, data_type, label_suffix,
                                 path_field = "path", enrich = NULL) {
  ctrl_insts  <- sample_data_instances(control_row, data_type)
  treat_insts <- sample_data_instances(treatment_row, data_type)
  ctrl_names  <- vapply(ctrl_insts, `[[`, character(1), "name")
  treat_names <- vapply(treat_insts, `[[`, character(1), "name")
  common <- intersect(ctrl_names, treat_names)
  if (length(common) == 0) return(list())

  safe_c <- gsub("[^A-Za-z0-9._-]", "_", control_row$display)
  safe_t <- gsub("[^A-Za-z0-9._-]", "_", treatment_row$display)
  ctrl_field  <- paste0("control_", path_field)
  treat_field <- paste0("treatment_", path_field)

  items <- lapply(common, function(nm) {
    safe_n     <- gsub("[^A-Za-z0-9._-]", "_", nm)
    ctrl_path  <- ctrl_insts[[which(ctrl_names == nm)[1]]]$path
    treat_path <- treat_insts[[which(treat_names == nm)[1]]]$path
    item <- list(
      key = paste0(viz, "_", safe_n, "_", safe_c, "_vs_", safe_t),
      label = paste0(nm, label_suffix, treatment_row$display, " (Cmp) vs ", control_row$display, " (Ref)"),
      type = viz, sample = NA_character_, file = NA_character_, colour_key = NA_character_,
      name = nm, samples = c(control_row$display, treatment_row$display),
      control_label = control_row$display, treatment_label = treatment_row$display
    )
    item[[ctrl_field]]  <- ctrl_path
    item[[treat_field]] <- treat_path
    item
  })

  if (!is.null(enrich)) items <- enrich(items, control_row, treatment_row)
  items
}

enrich_junctions_tpm_scaled <- function(items, control_row, treatment_row) {
  ctrl_quant  <- sample_data_path(control_row, "quantdir")
  treat_quant <- sample_data_path(treatment_row, "quantdir")
  if (is.na(ctrl_quant) || is.na(treat_quant)) return(list())
  lapply(items, function(it) modifyList(it, list(control_quant = ctrl_quant, treatment_quant = treat_quant)))
}

FIXED_ITEMS <- list(
  genes                 = "Gene annotation",
  protein_domains       = "Canonical protein annotation",
  transcripts_all       = "Transcripts",
  transcripts           = "Transcripts (filtered)",
  transcripts_annotated = "Transcripts annotated (filtered)"
)

fixed_item <- function(viz) {
  list(key = viz, label = FIXED_ITEMS[[viz]], type = viz,
       sample = NA_character_, file = NA_character_, colour_key = NA_character_)
}

detect_tpm_usage <- function(sample_row) {
  junc  <- sample_data_path(sample_row, "paired-end", name = "Junctions")
  quant <- sample_data_path(sample_row, "quantdir")
  if (is.null(quant) || is.na(quant)) return(NULL)
  safe_s <- gsub("[^A-Za-z0-9._-]", "_", sample_row$display)
  list(
    key = paste0("tpm_usage_", safe_s),
    label = paste0("Transcript counts · ", sample_row$display),
    type = "tpm_usage", sample = sample_row$display, file = NA_character_, colour_key = sample_row$display,
    junction_path = junc, quant_path = quant,
    broken = !file.exists(quant)
  )
}

build_track_items <- function(samples_list) {
  items <- list()
  for (row in samples_list) {
    for (item in detect_generic(row, "bw", "bw", path_field = "file")) items[[item$key]] <- item
    for (item in detect_generic(row, "hic", "hic", path_field = "file")) items[[item$key]] <- item
    for (item in detect_generic(row, "bed", "bed", path_field = "file")) items[[item$key]] <- item
    for (item in detect_generic(row, "sashimi", "paired-end", path_field = "junction_path")) items[[item$key]] <- item
    item <- detect_tpm_usage(row)
    if (!is.null(item)) items[[item$key]] <- item
  }
  for (viz in names(FIXED_ITEMS)) {
    item <- fixed_item(viz)
    items[[item$key]] <- item
  }
  items
}

build_comparison_track_items <- function(samples_list) {
  items <- list()
  n <- length(samples_list)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      control   <- samples_list[[i]]
      treatment <- samples_list[[j]]
      if (identical(control$display, treatment$display)) next
      for (item in detect_pair_generic(control, treatment, "bw_diff_stacked", "bw", " diff (stacked) · ", path_field = "path")) items[[item$key]] <- item
      for (item in detect_pair_generic(control, treatment, "junctions_mirrored_diff", "paired-end", " mirrored diff · ", path_field = "path")) items[[item$key]] <- item
      for (item in detect_pair_generic(control, treatment, "junctions_mirrored_diff_tpm_scaled", "paired-end", " mirrored diff (TPM-scaled) · ", path_field = "path", enrich = enrich_junctions_tpm_scaled)) items[[item$key]] <- item
    }
  }
  items
}

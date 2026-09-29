load_hic <- function(entry, region_params, settings) {
  is_cool <- tolower(tools::file_ext(entry$path)) %in% c("cool", "mcool")

  avail <- sort(if (is_cool) plotgardener::readCoolBpResolutions(entry$path) else strawr::readHicBpResolutions(entry$path))
  res <- if (identical(entry$hic_res_mode %||% settings$hic_res_mode %||% "auto", "finest")) {
    min(avail)
  } else {
    cands <- avail[((region_params$end - region_params$start) / avail) >= (settings$hic_max_bins_per_view %||% 500)]
    if (length(cands) > 0) max(cands) else avail[1]
  }

  # hic's chrom vocab is keyed BY RESOLUTION (a coarser .mcool can drop small
  # contigs), so this can't be resolved before `res` is picked above.
  vocab <- if (is_cool) {
    plotgardener::readCoolChroms(entry$path, resolution = res)$name
  } else {
    v <- strawr::readHicChroms(entry$path)$name
    v[v != "All"]
  }
  chrom <- resolve_chrom(region_params$chrom, vocab)
  if (is.null(chrom)) {
    entry$data <- list(matrix = NULL, chrom = NA_character_, resolution = res, failed = TRUE, no_chrom = TRUE)
    return(entry)
  }

  norm <- settings$hic_norm %||% "KR"
  # format(scientific=FALSE): paste0() on a large round double (e.g. 80000000)
  # silently renders as "8e+07", which strawr::straw() parses as an empty range.
  locus <- paste0(chrom, ":", format(region_params$start, scientific = FALSE), ":", format(region_params$end, scientific = FALSE))
  m <- if (is_cool) {
    plotgardener::readCool(file = entry$path, chrom = chrom, chromstart = region_params$start, chromend = region_params$end, resolution = res, norm = norm)
  } else {
    strawr::straw(norm, entry$path, locus, locus, "BP", res)
  }
  colnames(m) <- c("x", "y", "counts")

  entry$data <- list(matrix = m, chrom = chrom, resolution = res)
  entry
}

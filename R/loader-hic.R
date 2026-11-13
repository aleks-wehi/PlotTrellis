load_hic <- function(entry, region_params, settings) {
  .pt_check_file(entry$path)
  is_cool <- tolower(tools::file_ext(entry$path)) %in% c("cool", "mcool")

  avail <- tryCatch(if (is_cool) plotgardener::readCoolBpResolutions(entry$path) else strawr::readHicBpResolutions(entry$path),
                    error = function(e) .pt_read_error(entry$path, e))
  avail <- sort(avail)
  res <- if (identical(settings$hic_res_mode, "finest")) {
    min(avail)
  } else if (identical(settings$hic_res_mode, "fixed")) {
    # Exactly the bin size asked for, or an error naming the closest ones: a different size
    # would silently draw something other than what was asked for. The track error handler
    # in render_all() turns the error into this track's in-figure message.
    want <- clean_optional_numeric(settings$hic_resolution)
    if (is.null(want)) stop("Hi-C resolution mode is Fixed but no bin size is set.", call. = FALSE)
    if (!(want %in% avail)) {
      near <- avail[order(abs(avail - want), avail)][seq_len(min(3, length(avail)))]
      stop("this file has no ", format(want, big.mark = ",", scientific = FALSE), " bp bins; closest: ",
           paste(format(sort(near), big.mark = ",", scientific = FALSE, trim = TRUE), collapse = ", "), " bp.", call. = FALSE)
    }
    want
  } else {
    cands <- avail[((region_params$end - region_params$start) / avail) >= settings$hic_max_bins_per_view]
    if (length(cands) > 0) max(cands) else avail[1]
  }

  # hic's chrom vocab is keyed BY RESOLUTION (a coarser .mcool can drop small
  # contigs), so this can't be resolved before `res` is picked above.
  vocab <- tryCatch(if (is_cool) {
    plotgardener::readCoolChroms(entry$path, resolution = res)$name
  } else {
    v <- strawr::readHicChroms(entry$path)$name
    v[v != "All"]
  }, error = function(e) .pt_read_error(entry$path, e))
  chrom <- resolve_chrom(region_params$chrom, vocab)
  if (is.null(chrom)) {
    entry$data <- list(matrix = NULL, chrom = NA_character_, resolution = res, failed = TRUE, no_chrom = TRUE)
    return(entry)
  }

  norm <- settings$hic_norm
  # format(scientific=FALSE): paste0() on a large round double (e.g. 80000000)
  # silently renders as "8e+07", which strawr::straw() parses as an empty range.
  locus <- paste0(chrom, ":", format(region_params$start, scientific = FALSE), ":", format(region_params$end, scientific = FALSE))
  # A .hic that lists the chromosome but holds no contacts on it (chrM, say) has nothing in view:
  # an empty matrix, which render_hic() draws from instead of the file (`no_contacts`).
  no_contacts <- FALSE
  m <- tryCatch(if (is_cool) {
    plotgardener::readCool(file = entry$path, chrom = chrom, chromstart = region_params$start, chromend = region_params$end, resolution = res, norm = norm)
  } else {
    strawr::straw(norm, entry$path, locus, locus, "BP", res)
  }, error = function(e) {
    if (!grepl("doesn't have the given chr_chr map", conditionMessage(e), fixed = TRUE)) .pt_read_error(entry$path, e)
    no_contacts <<- TRUE
    data.frame(x = integer(0), y = integer(0), counts = numeric(0))
  })
  colnames(m) <- c("x", "y", "counts")

  entry$data <- list(matrix = m, chrom = chrom, resolution = res, no_contacts = no_contacts)
  entry
}

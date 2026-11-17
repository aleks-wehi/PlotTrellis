# Normalisation names offered in the settings: what Juicer writes into a .hic
# (KR, SCALE, VC, VC_SQRT), cooler's own balancing in a .cool/.mcool (BALANCE,
# its "weight" column) and NONE (raw counts, in every file). Others can be typed
# in.
HIC_NORMS <- c("KR", "SCALE", "VC", "VC_SQRT", "BALANCE", "NONE")

# plotgardener draws at most this many bins across a view; for more it re-reads
# the file at a bin size of its own (hic_limit()), which our colour range would
# not match.
HIC_MAX_BINS <- 1000

# The bin size a Hi-C track is read and drawn at, from the file's bin sizes
# `avail` (sorted) and the view width in bp. Only sizes giving at most
# HIC_MAX_BINS bins across the view are used (the coarsest size if none does),
# so plotgardener draws at the size the colour range was taken from.
.pt_hic_pick_resolution <- function(avail, width, settings) {
  fits <- avail[width / avail <= HIC_MAX_BINS]
  if (!length(fits)) fits <- max(avail)
  if (identical(settings$hic_res_mode, "finest")) return(min(fits))
  if (identical(settings$hic_res_mode, "fixed")) {
    # Exactly the bin size asked for, or an error naming the closest ones: a
    # different size would silently draw something other than what was asked
    # for. The track error handler in render_all() turns the error into this
    # track's in-figure message.
    want <- clean_optional_numeric(settings$hic_resolution)
    if (is.null(want))
      stop("Hi-C Resolution is Fixed but no Bin size (bp) is set.",
           call. = FALSE)
    if (!(want %in% avail)) {
      near <- avail[order(abs(avail - want), avail)][
        seq_len(min(3, length(avail)))
      ]
      stop("this file has no ",
           format(want, big.mark = ",", scientific = FALSE),
           " bp bins; closest: ",
           paste(format(sort(near), big.mark = ",", scientific = FALSE,
                        trim = TRUE), collapse = ", "),
           " bp.", call. = FALSE)
    }
    # Too fine to draw across this view: the finest size that fits, as
    # plotgardener would draw a coarser one anyway. The bin size drawn is shown
    # under the colour key.
    return(if (want %in% fits) want else min(fits))
  }
  # Auto: the coarsest size that still gives at least `hic_max_bins_per_view`
  # bins (the finest that fits if none does).
  min_bins <- clean_optional_numeric(settings$hic_max_bins_per_view) %||% 500
  cands <- fits[width / fits >= min_bins]
  if (length(cands)) max(cands) else min(fits)
}

load_hic <- function(entry, region_params, settings) {
  .pt_check_file(entry$path)
  is_cool <- tolower(tools::file_ext(entry$path)) %in% c("cool", "mcool")

  avail <- tryCatch(
    if (is_cool) plotgardener::readCoolBpResolutions(entry$path) else
      strawr::readHicBpResolutions(entry$path),
    error = function(e) .pt_read_error(entry$path, e)
  )
  avail <- sort(avail)
  res <- .pt_hic_pick_resolution(avail,
                                 region_params$end - region_params$start,
                                 settings)

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
    entry$data <- list(matrix = NULL, chrom = NA_character_, resolution = res,
                       failed = TRUE, no_chrom = TRUE)
    return(entry)
  }

  norm <- settings$hic_norm
  # format(scientific=FALSE): paste0() on a large round double (e.g. 80000000)
  # silently renders as "8e+07", which strawr::straw() parses as an empty range.
  locus <- paste0(chrom, ":", format(region_params$start, scientific = FALSE),
                  ":", format(region_params$end, scientific = FALSE))
  # A .hic that lists the chromosome but holds no contacts on it (chrM, say) has
  # nothing in view: an empty matrix, which render_hic() draws from instead of
  # the file (`no_contacts`).
  no_contacts <- FALSE
  m <- tryCatch(if (is_cool) {
    plotgardener::readCool(file = entry$path, chrom = chrom,
                           chromstart = region_params$start,
                           chromend = region_params$end, resolution = res,
                           norm = norm)
  } else {
    strawr::straw(norm, entry$path, locus, locus, "BP", res)
  }, error = function(e) {
    # A .cool contig with no contacts at the end of its chromosome list:
    # readCool() runs past the end of the file's index. Said in the track
    # (render_all()'s error handler), not a read error.
    if (is_cool &&
        grepl("index exceeds HDF5-array dimension", conditionMessage(e),
              fixed = TRUE))
      stop("this file has no contacts on ", chrom, ".", call. = FALSE)
    # The Normalisation set is not in this file (at this bin size, for a .hic):
    # say so and what to pick, instead of the reader's own message, which is
    # long and gets cut.
    if (grepl("is not found in|did not contain .*normali[sz]ation",
              conditionMessage(e))) {
      has <- tryCatch(
        unique(unlist(if (is_cool) plotgardener::readCoolNorms(entry$path) else
                        strawr::readHicNormTypes(entry$path))),
        error = function(e2) character(0)
      )
      stop("this file has no ", norm, " normalisation",
           if (!is_cool) paste0(" at ",
                                format(res, big.mark = ",", scientific = FALSE),
                                " bp bins"),
           if (length(has)) paste0("; it has ", paste(has, collapse = ", "),
                                   ": pick one of those.") else ": try NONE.",
           call. = FALSE)
    }
    if (!grepl("doesn't have the given chr_chr map", conditionMessage(e),
               fixed = TRUE))
      .pt_read_error(entry$path, e)
    no_contacts <<- TRUE
    data.frame(x = integer(0), y = integer(0), counts = numeric(0))
  })
  colnames(m) <- c("x", "y", "counts")
  # Contacts in view but none with a value: the normalisation has no vectors at
  # this bin size (a sparse file's finest bins, say), which would draw an empty
  # panel with no explanation. The track error handler in render_all() shows
  # this as the track's message instead. Only when this bin size is the one
  # drawn: over HIC_MAX_BINS across the view (only when even the file's coarsest
  # size is too fine), plotgardener's hic_limit() re-reads the file at another
  # size.
  drawn_at_res <- (region_params$end - region_params$start) / res <=
    HIC_MAX_BINS
  if (!no_contacts && drawn_at_res && nrow(m) > 0 &&
      !any(is.finite(m$counts))) {
    bp <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    normalised <- !identical(toupper(norm), "NONE")
    coarser <- avail[avail > res]
    # Kept under the 140 characters render_all() shows without cutting.
    fixes <- c(if (normalised) "set Hi-C normalisation to NONE",
               if (length(coarser)) paste0("Resolution Fixed with a ",
                                           bp(coarser[1]), " bp bin size"))
    stop("no ", if (normalised) paste(norm, "values") else "contacts", " at ",
         bp(res), " bp bins in this view",
         if (length(fixes)) paste0(": ", paste(fixes, collapse = ", or ")),
         ".", call. = FALSE)
  }

  entry$data <- list(matrix = m, chrom = chrom, resolution = res,
                     no_contacts = no_contacts)
  entry
}

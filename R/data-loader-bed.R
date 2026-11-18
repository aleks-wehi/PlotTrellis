# Reads the file once (.pt_read_bed()) for its chromosome names; render_bed()
# draws the same rows. A file that is missing or can't be read stops; an empty
# file, or one with only header, comment and blank lines, is a track with
# nothing in view (`empty`).
load_bed <- function(entry, region_params) {
  .pt_check_file(entry$path)
  bed <- .pt_read_bed(entry$path)
  if (nrow(bed) == 0) {
    entry$data <- list(no_chrom = FALSE, empty = TRUE)
    return(entry)
  }
  chrom <- resolve_chrom(region_params$chrom, unique(bed[[1]]))
  if (is.null(chrom)) {
    entry$data <- list(no_chrom = TRUE)
    return(entry)
  }
  entry$data <- list(no_chrom = FALSE, chrom = chrom, bed = bed)
  entry
}

# A BED file's intervals: chromosome, start and end (V1-V3; plotRanges() draws
# nothing else of them). Dropped: blank lines, header and comment lines (`#...`,
# `track`, `browser`), and the lines before the first one with a numeric start
# and end (a header row without '#', e.g. "chrom start end"). No line is lost to
# fread()'s guesses about the file's layout: it reads the file at once (a
# compressed file too) from the first interval, and its result is kept only if
# it read every row cleanly; a file with a `track`, `#` or tab-only line further
# down, or spaces for tabs, is read line by line instead. A file whose lines
# have fewer than 3 fields, or a start or end that isn't a number, stops with
# the reason; a single header row and nothing else, or only blank and comment
# lines, is a file with no intervals.
.pt_read_bed <- function(path) {
  none <- data.frame(V1 = character(0), V2 = numeric(0), V3 = numeric(0))
  if (file.size(path) == 0) return(none)
  fields <- function(x, k) {
    has <- grepl(paste0("^\\s*", strrep("\\S+\\s+", k - 1), "\\S+"), x,
                 perl = TRUE)
    out <- rep(NA_character_, length(x))
    out[has] <- sub(paste0("^\\s*", strrep("\\S+\\s+", k - 1), "(\\S+).*$"),
                    "\\1", x[has], perl = TRUE)
    out
  }
  number <- function(x) suppressWarnings(as.numeric(x))
  interval <- function(x) !grepl("^\\s*(#|(track|browser)(\\s|$))", x) &
    !is.na(number(fields(x, 2))) &
    !is.na(number(fields(x, 3)))
  # The fast way: fread() from the first interval among the first lines.
  # (A file that can't be opened here gets fread()'s own message below.)
  top <- tryCatch(suppressWarnings({
    con <- gzfile(path, "r")
    tryCatch(readLines(con, n = 100L, warn = FALSE), finally = close(con))
  }), error = function(e) character(0))
  # a byte order mark, which fread() drops too
  top <- sub("^\ufeff", "", top, useBytes = TRUE)
  first <- match(TRUE, interval(top))
  if (!is.na(first)) {
    lead <- c(fields(top[first], 1), fields(top[first], 2),
              fields(top[first], 3))
    warned <- FALSE
    d <- tryCatch(withCallingHandlers(
      .pt_fread(path,
                sep = if (grepl("\t", top[first], fixed = TRUE)) "\t" else " ",
                header = FALSE, skip = first - 1L, fill = TRUE,
                blank.lines.skip = TRUE, quote = "",
                colClasses = list(character = 1L)),
      warning = function(w) {
        warned <<- TRUE; invokeRestart("muffleWarning")
      }), error = function(e) NULL)
    # Kept only if read cleanly, and from that same first line (fread() counts
    # lines as readLines() does except for odd line ends, such as "\r\r\n").
    if (!warned && !is.null(d) && nrow(d) > 0 && ncol(d) >= 3 &&
        is.numeric(d[[2]]) && is.numeric(d[[3]]) &&
        !anyNA(d[[2]]) && !anyNA(d[[3]]) &&
        identical(as.character(d[[1]][1]), lead[1]) &&
        isTRUE(d[[2]][1] == number(lead[2])) &&
        isTRUE(d[[3]][1] == number(lead[3]))) {
      d <- as.data.frame(d[, 1:3])
      # A commented-out interval ("#chr1 ...") or a track line with numbers
      # further down.
      v1 <- as.character(d[[1]])
      out <- startsWith(v1, "#") | v1 %in% c("track", "browser")
      return(if (any(out)) d[!out, , drop = FALSE] else d)
    }
  }
  # Line by line.
  l <- tryCatch(.pt_fread(path, sep = "\n", header = FALSE, quote = "",
                          colClasses = "character")[[1]],
                error = function(e) if (grepl("Input is either empty",
                                              conditionMessage(e),
                                              fixed = TRUE)) character(0)
                                    else .pt_read_error(path, e))
  ln <- seq_along(l)  # file line numbers, for messages
  i <- which(substr(l, 1, 1) %in% c("#", "t", "b", "", " ", "\t"))
  drop <- i[grepl("^\\s*(#|(track|browser)(\\s|$)|$)", l[i])]
  if (length(drop)) { l <- l[-drop]; ln <- ln[-drop] }
  if (!length(l)) return(none)
  # The first three fields of each line that has them, split on spaces and tabs.
  three <- grepl("^\\s*\\S+\\s+\\S+\\s+\\S+", l, perl = TRUE)
  f <- list(rep(NA_character_, length(l)), rep(NA_character_, length(l)),
            rep(NA_character_, length(l)))
  if (any(three)) {
    parts <- strsplit(sub("^\\s*(\\S+)\\s+(\\S+)\\s+(\\S+).*$", "\\1\t\\2\t\\3",
                          l[three], perl = TRUE),
                      "\t", fixed = TRUE)
    # exactly 3 fields each
    parts <- matrix(unlist(parts, use.names = FALSE), nrow = 3L)
    for (k in 1:3) f[[k]][three] <- parts[k, ]
  }
  start <- number(f[[2]])
  end <- number(f[[3]])
  ok <- !is.na(start) & !is.na(end)
  first <- match(TRUE, ok)
  n_fields <- function(x) lengths(strsplit(trimws(x), "\\s+"))
  if (is.na(first)) {
    n <- n_fields(l)
    if (length(l) == 1 && n >= 3) return(none)  # a header row alone
    if (all(n < 3))
      stop(basename(path), " is not a BED file: it has ", max(n), " column",
           if (max(n) != 1) "s",
           ", and needs at least 3 (chromosome, start, end)", call. = FALSE)
    eg <- min(2L, length(l))  # past a header row
    stop(basename(path),
         " is not a BED file: no line has numbers for start and end (columns 2 and 3), ",
         "e.g. line ", ln[eg], ": ", substr(l[eg], 1, 60), call. = FALSE)
  }
  keep <- seq.int(first, length(l))
  bad <- keep[match(FALSE, ok[keep])]
  if (!is.na(bad)) {
    n <- n_fields(l[bad])
    stop(basename(path), " is not a BED file: line ", ln[bad], " ",
         if (n < 3) paste0("has ", n, " field", if (n != 1) "s",
                           ", and needs at least 3 (chromosome, start, end)")
         else "has a start or end that isn't a number",
         ": ", substr(l[bad], 1, 60), call. = FALSE)
  }
  data.frame(V1 = f[[1]][keep], V2 = start[keep], V3 = end[keep])
}

# data.table::fread(), for a gzip or bzip2 file too. fread() decompresses those
# only through the R.utils package, which PlotTrellis does not import, so this
# decompresses them with base R into a temporary file first, chosen as fread()
# chooses: a .gz/.bgz name or gzip's first bytes with gzfile(), else a .bz2 name
# or bzip2's with bzfile() (both also read a file that turns out to be plain
# text).
.pt_fread <- function(path, ...) {
  magic <- readBin(path, "raw", 3L)
  gz <- grepl("\\.b?gz$", path) || identical(magic[1:2], as.raw(c(0x1f, 0x8b)))
  bz <- !gz && (grepl("\\.bz2$", path) || identical(magic, charToRaw("BZh")))
  if (!gz && !bz) return(data.table::fread(path, ...))
  tmp <- tempfile()
  on.exit(unlink(tmp))
  src <- (if (gz) gzfile else bzfile)(path, "rb")
  dst <- file(tmp, "wb")
  tryCatch(repeat {
    chunk <- readBin(src, "raw", 1e7)
    if (length(chunk) == 0) break
    writeBin(chunk, dst)
  }, finally = { close(src); close(dst) })
  data.table::fread(tmp, ...)
}

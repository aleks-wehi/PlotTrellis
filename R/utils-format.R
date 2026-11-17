# Number, coordinate and label formatting for the composer (R/app.R), as plain
# functions. .pt_fmt_num() rounds before formatting, unlike .fmt_num() in
# R/boilerplate.R, which the figures use: the two are not interchangeable.

# Bolds whatever precedes the first " · " in a label (or the whole label if it
# has none); format-driven, not an enumerated word list, since a track's label
# can start with any user-given instance name.
.pt_format_track_label <- function(label) {
  sub("^([^\u00B7]+?)( \u00B7 |$)", "<strong>\\1</strong>\\2", label)
}

# Main chromosomes first (chr1..chrN, chrX, chrY, chrM), then every
# alt/random/fix/chrUn contig alphabetically below them -- kept, not filtered.
.pt_sort_chrs <- function(chrs) {
  is_std <- base::grepl("^chr[0-9]+$", chrs)
  nums   <- suppressWarnings(as.integer(base::sub("^chr", "", chrs)))
  sex_mt <- intersect(c("chrX", "chrY", "chrM"), chrs)
  c(chrs[is_std][order(nums[is_std])], sex_mt,
    sort(setdiff(chrs[!is_std], sex_mt)))
}
.pt_fmt_num <- function(x, digits = 2) {
  s <- formatC(round(x, digits), format = "f", digits = digits)
  s <- sub("0+$", "", s); sub("\\.$", "", s)
}
.pt_format_bp <- function(bp) {
  if (bp >= 1e6) paste0(.pt_fmt_num(bp / 1e6), "Mb")
  else if (bp >= 1e3) paste0(.pt_fmt_num(bp / 1e3), "Kb")
  else paste0(bp, "bp")
}
.pt_format_coord <- function(n) {
  if (n >= 1e6) paste0(.pt_fmt_num(n / 1e6, 3), "Mb")
  else if (n >= 1e3) paste0(.pt_fmt_num(n / 1e3, 3), "Kb")
  else paste0(n, " bp")
}

# Parse human-readable bp strings: "17000000", "17Mb", "17.56Mb", "15kb",
# "150KB"
.pt_parse_bp <- function(s) {
  s <- trimws(s)
  if (grepl("^[0-9]+$", s)) return(as.integer(s))
  m <- regmatches(s, regexec("^([0-9]*\\.?[0-9]+)\\s*(kb|mb)$", s,
                             ignore.case = TRUE))[[1]]
  if (length(m) < 3) return(NA_integer_)
  num  <- as.numeric(m[2])
  unit <- toupper(m[3])
  if (unit == "KB") as.integer(round(num * 1e3))
  else if (unit == "MB") as.integer(round(num * 1e6))
  else NA_integer_
}

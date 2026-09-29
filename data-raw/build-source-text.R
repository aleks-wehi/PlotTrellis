# Regenerate R/sysdata.rda's PT_SOURCE_TEXT -- run this after adding or
# editing any function in R/*.R. Captures each top-level function's exact
# original source text (comments, formatting) here, at build time, because
# an installed package doesn't retain function srcrefs by default (needs
# R_KEEP_PKG_SOURCE=yes/--with-keep.source, which no ordinary install sets)
# -- so assemble_download_script() can't recover it from the live object at
# runtime the way it could when the package is loaded from source. See
# .pt_source_text() in R/download-script.R, the only reader of this data.

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)

# Each function-defining statement is eval()'d alone (never the whole file --
# some R/ files build top-level constants via Imports this dev environment
# may not have installed) into a shared env, so later definitions can
# reference earlier ones by name inside their own bodies without those
# bodies ever actually running here. Evaluating just `name <- function(...)
# ...` attaches the srcref of the function literal itself (comments,
# formatting, but NOT the "name <-" part) to the resulting closure -- the
# same attribute .pt_source_text() used to read straight off the live,
# installed-with-keep.source object.
env <- new.env()

PT_SOURCE_TEXT <- list()
for (f in r_files) {
  exprs <- parse(f, keep.source = TRUE)
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (!is.call(e) || length(e) < 3) next
    op <- as.character(e[[1]])
    if (!op %in% c("<-", "=")) next
    target <- e[[2]]
    rhs <- e[[3]]
    if (!is.symbol(target) || !is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
    eval(e, envir = env)
    fn <- get(as.character(target), envir = env)
    PT_SOURCE_TEXT[[as.character(target)]] <- paste(as.character(attr(fn, "srcref")), collapse = "\n")
  }
}

save(PT_SOURCE_TEXT, file = "R/sysdata.rda")

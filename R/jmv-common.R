# jmv-common.R -- helpers shared by the jamovi analyses (the *.b.R files).
#
# The jamovi layer is deliberately thin: it turns options into arguments of
# the engine functions (invarianceWiz(), cfaWiz()), caches what they return
# and copies the returned data frames into jamovi tables. No statistics here.

# The `factors` Array option -> named list of item vectors (empty boxes dropped).
psw_jmv_factors <- function(factors) {
  out <- list()
  for (k in seq_along(factors)) {
    vars <- unlist(factors[[k]]$vars)
    if (!length(vars)) next
    label <- factors[[k]]$label
    if (is.null(label) || !nzchar(trimws(label))) label <- paste("Factor", k)
    out[[length(out) + 1L]] <- vars
    names(out)[length(out)] <- trimws(label)
  }
  if (length(out)) names(out) <- make.unique(names(out), sep = " ")
  out
}

# Enough has been specified to try fitting (the engine does the real checks).
psw_jmv_ready <- function(factors) {
  length(factors) > 0L && all(lengths(factors) >= 2L) && sum(lengths(factors)) >= 3L
}

# jamovi shows NA as "NaN"; an empty string gives an empty cell instead.
psw_jmv_val <- function(x) {
  if (length(x) != 1L || is.na(x)) '' else x
}

# Copy rows (a list of named lists) into a table that starts with zero rows.
# Works both on a fresh table and on one jamovi restored from a previous run.
psw_jmv_fill <- function(table, rows) {
  for (i in seq_along(rows)) {
    values <- lapply(rows[[i]], psw_jmv_val)
    if (i <= table$rowCount) table$setRow(rowNo = i, values = values)
    else table$addRow(rowKey = i, values = values)
  }
}

psw_jmv_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# The message box above the tables: instructions, the plan, or nothing.
psw_jmv_todo <- function(html, lines = NULL, notes = NULL) {
  if (is.null(lines)) {
    html$setVisible(FALSE)
    return(invisible())
  }
  body <- paste0("<p>", lines, "</p>", collapse = "")
  if (length(notes))
    body <- paste0(body, "<ul>", paste0("<li>", psw_jmv_escape(notes), "</li>", collapse = ""), "</ul>")
  html$setContent(body)
  html$setVisible(TRUE)
}

# Engine notes (data issues, lavaan warnings) become footnotes of a table.
psw_jmv_notes <- function(table, notes, prefix = "note") {
  for (i in seq_along(notes)) table$setNote(paste0(prefix, i), notes[i], init = FALSE)
}

# "scaled chi-square; robust CFI, TLI and RMSEA" for a table footnote.
psw_jmv_variant_text <- function(variant) {
  if (all(is.na(variant))) return(NULL)
  chi <- switch(variant[["chisq"]], scaled = "\u03C7\u00B2 is the scaled (robust) test statistic", "\u03C7\u00B2 is the standard test statistic")
  idx <- if (identical(variant[["cfi"]], variant[["rmsea"]]))
    sprintf("CFI, TLI and RMSEA are the %s versions", variant[["cfi"]])
  else sprintf("CFI and TLI are the %s versions, RMSEA the %s version", variant[["cfi"]], variant[["rmsea"]])
  paste0(chi, "; ", idx, ".")
}

psw_jmv_time_text <- function(info) {
  sprintf("lavaan %s; fitted in %.1f s%s.", info$lavaan, info$time,
          if (info$workers > 0L) sprintf(" on %d processor cores", info$workers) else "")
}

# ---- figures ----------------------------------------------------------------
# .run() stores what the figure needs as the image's state and sets its size;
# the render function of every analysis (.diagram) is psw_jmv_diagram().
psw_jmv_set_figure <- function(image, loadings, cors, panel_col, type, highlight, caption = NULL) {
  size <- if (type == "comb") psw_comb_size(loadings, cors, panel_col)
          else psw_diagram_size(loadings, cors, panel_col, tiles = type == "tiles")
  image$setSize(size[["width"]], size[["height"]])
  image$setState(list(loadings = loadings, cors = cors, panel_col = panel_col, type = type,
                      highlight = isTRUE(highlight), caption = caption))
}

psw_jmv_diagram <- function(image) {
  st <- image$state
  if (is.null(st)) return(FALSE)
  print(psw_figure(st$loadings, st$cors, st$panel_col, st$type, st$highlight, st$caption))
  TRUE
}

# Criteria and cutoffs from the shared "Decision Criteria" options.
psw_jmv_criteria <- function(opt) {
  list(criteria = c("cfi", "rmsea", "srmr", "chisq")[c(opt$useCfi, opt$useRmsea, opt$useSrmr, opt$useChisq)],
       cutoffs = list(cfi = opt$cutCfi, rmsea = opt$cutRmsea, srmr.loadings = opt$cutSrmrLoad,
                      srmr = opt$cutSrmr, alpha = opt$alpha))
}

# Rows of a model-fit table from an engine `fit` data frame.
psw_jmv_fit_rows <- function(f) {
  lapply(seq_len(nrow(f)), function(i)
    list(model = f$label[i], chisq = f$chisq[i], df = round(f$df[i]), p = f$pvalue[i],
         cfi = f$cfi[i], tli = f$tli[i], rmsea = f$rmsea[i],
         rmseaLo = f$rmsea.lo[i], rmseaHi = f$rmsea.hi[i], srmr = f$srmr[i]))
}

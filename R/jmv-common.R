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
# .run() stores WHAT to draw as the image's state (a kind and the arguments of
# its drawing function, never a plot object: the state is saved with the
# analysis) and sets the image size the layout expects. Every analysis uses
# psw_jmv_diagram() as its render function.
psw_jmv_build <- function(kind, args) {
  do.call(switch(kind, figure = psw_figure, scalogram = psw_scalogram_plot, antenna = psw_antenna,
                 facing = psw_facing, bfpath = psw_bifactor_path, brush = psw_brush,
                 means = psw_plot_means), args)
}

psw_jmv_plot <- function(image, kind, args) {
  size <- attr(psw_jmv_build(kind, args), "psw_size")
  image$setSize(size[["width"]], size[["height"]])
  image$setState(list(kind = kind, args = args))
}

psw_jmv_diagram <- function(image) {
  st <- image$state
  if (is.null(st)) return(FALSE)
  print(psw_jmv_build(st$kind, st$args))
  TRUE
}

# The comb plot and the path diagram of a simple-structure model, from the
# shared figure options (showComb, combThick, showPath, pathNumbers, pathWidth, plotFit, showPlotGuide).
psw_jmv_figures <- function(out, opt, loadings, cors, panel_col = "group", caption = NULL, resid = NULL, higher = NULL,
                            fit = NULL, estimator = NULL) {
  fit_text <- if (!is.null(opt$plotFit) && opt$plotFit != "none" && !is.null(fit)) psw_format_fit_line(fit, estimator) else NULL
  fit_pos <- if (!is.null(opt$plotFit)) opt$plotFit else "none"
  show_guide <- !isFALSE(opt$showPlotGuide)
  if (isTRUE(opt$showComb)) {
    psw_jmv_plot(out$comb, "figure", list(loadings = loadings, cors = cors, panel_col = panel_col, type = "comb",
                                          caption = caption, opts = list(thick = opt$combThick, resid = resid, higher = higher,
                                                                         fit_text = fit_text, fit_pos = fit_pos)))
    if (!is.null(out$combDesc) && show_guide)
      out$combDesc$setContent(psw_jmv_comb_desc(thick = isTRUE(opt$combThick),
                                                has_resid = !is.null(resid) && nrow(resid) > 0L,
                                                has_higher = !is.null(higher),
                                                has_cors = !is.null(cors) && nrow(cors) > 0L))
  }
  if (isTRUE(opt$showPath))
    psw_jmv_plot(out$pathFig, "figure", list(loadings = loadings, cors = cors, panel_col = panel_col, type = "path",
                                          caption = caption, opts = list(numbers = opt$pathNumbers, width = opt$pathWidth,
                                                                         fit_text = fit_text, fit_pos = fit_pos)))
}

# Flag poor model-fit indices in red when the user turns on "Flag weak / misfitting values".
psw_jmv_flag_fit <- function(table, f) {
  if (is.null(f) || !nrow(f)) return(invisible())
  for (i in seq_len(nrow(f))) {
    if ("cfi" %in% names(f) && isTRUE(f$cfi[i] < 0.90))
      table$addFormat(rowNo = i, col = "cfi", format = jmvcore::Cell.NEGATIVE)
    if ("tli" %in% names(f) && isTRUE(f$tli[i] < 0.90))
      table$addFormat(rowNo = i, col = "tli", format = jmvcore::Cell.NEGATIVE)
    if ("rmsea" %in% names(f) && isTRUE(f$rmsea[i] > 0.08))
      table$addFormat(rowNo = i, col = "rmsea", format = jmvcore::Cell.NEGATIVE)
    if ("srmr" %in% names(f) && isTRUE(f$srmr[i] > 0.08))
      table$addFormat(rowNo = i, col = "srmr", format = jmvcore::Cell.NEGATIVE)
  }
  table$setNote("hl_fit", "In red: CFI or TLI below .90, RMSEA or SRMR above .08.", init = FALSE)
}

# Compact HTML plot guides displayed above novel figures (toggleable via `showPlotGuide`).
psw_jmv_guide_box <- function(title, body) {
  paste0(
    "<div style='font-size: 12px; color: #3f444a; line-height: 1.45; ",
    "border-left: 3px solid #8a94a0; padding: 6px 10px; margin: 4px 0 8px 0; background: #f8f9fa;'>",
    "<b>", title, ":</b> ", body,
    "</div>"
  )
}

psw_jmv_comb_desc <- function(thick = FALSE, has_resid = FALSE, has_higher = FALSE, has_cors = FALSE) {
  extras <- character()
  if (has_cors && !has_higher) extras <- c(extras, "Arcs link correlated factors (thicker = stronger).")
  if (has_higher) extras <- c(extras, "The top tooth on each factor shows how much of it the second-order factor explains.")
  if (has_resid) extras <- c(extras, "Brackets at the right mark freed residual correlations.")
  if (thick) extras <- c(extras, "Tooth thickness reflects each item's observed variance.")
  body <- paste(c(
    "Each vertical spine is a factor and each horizontal tooth is one item's full variance (0 to 1).",
    "The <b>dark bar</b> is the share of variance explained by the factor (&lambda;&sup2;, labeled with standardized &lambda;); the <b>pale remainder</b> is unique/error variance.",
    "Longer dark teeth mean stronger items.",
    extras
  ), collapse = " ")
  psw_jmv_guide_box("Reading the Comb Plot", body)
}

psw_jmv_antenna_desc <- function() {
  psw_jmv_guide_box(
    "Reading the Antenna Plot",
    paste(
      "Each row is one item split at a central mast (0).",
      "The <b>dark bar to the left</b> is the share of item variance explained by the general factor (&lambda;<sub>g</sub>&sup2;);",
      "the <b>grey bar to the right</b> is the share explained by its specific factor (&lambda;<sub>s</sub>&sup2;), labeled with standardized loadings.",
      "Items leaning left mostly reflect the general trait; items leaning right retain strong subscale-specific variance."
    )
  )
}

psw_jmv_facing_desc <- function() {
  psw_jmv_guide_box(
    "Reading Facing Combs",
    paste(
      "Each row is one item's full variance (0 to 1).",
      "General variance grows from the <b>left (dark)</b> and specific variance grows from the <b>right (grey)</b>, labeled with standardized loadings.",
      "The <b>pale gap in the middle</b> is unexplained unique variance &mdash; combs that nearly meet indicate well-explained items."
    )
  )
}

psw_jmv_brush_desc <- function() {
  psw_jmv_guide_box(
    "Reading the Brush Plot",
    paste(
      "Each vertical handle is a factor, with horizontal bristles showing squared standardized loadings (&lambda;&sup2;, labeled when |&lambda;| &ge; .10).",
      "<b>Dark bristles</b> are target loadings on an item's own factor; <b>grey bristles</b> are cross-loadings on other factors; the final column is unique variance.",
      "A clean factor looks like a toothbrush &mdash; strong bristles on its own items and a bare handle elsewhere."
    )
  )
}

psw_jmv_scalogram_desc <- function(binary = FALSE) {
  psw_jmv_guide_box(
    "Reading the Scalogram",
    paste(
      "Respondents are sorted vertically from lowest to highest total score against the items.",
      if (binary)
        "In a perfect Guttman scale everything left of the step line is dark (endorsed) and everything right of it is light."
      else
        "In a clean cumulative scale the dark area grows smoothly from left to right as the total score rises, without scattered islands."
    )
  )
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

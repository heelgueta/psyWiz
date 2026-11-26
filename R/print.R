# print.R -- console output for R users.
#
# Keeps the look of the original cfawiz(): compact fixed-width tables with one
# column per estimator (or per step) and a red-to-blue colour ramp that makes
# weak loadings and poor fit stand out. Colours are off by default; ask for
# them with `print(x, colored = TRUE)` (used only if the terminal supports them).

# the ramps are in R/colour.R; the console keeps the original red-to-blue one
psw_palette <- psw_palettes$redblue

psw_use_color <- function(colored = TRUE) {
  if (!isTRUE(colored) || nzchar(Sys.getenv("NO_COLOR"))) return(FALSE)
  if (requireNamespace("cli", quietly = TRUE)) return(cli::num_ansi_colors() >= 256L)
  isatty(stdout()) && !identical(Sys.getenv("TERM"), "dumb")
}

psw_paint <- function(text, index) {
  rgb <- grDevices::col2rgb(psw_palette[max(1L, min(index, length(psw_palette)))])
  sprintf("\033[38;2;%d;%d;%dm%s\033[39m", rgb[1L], rgb[2L], rgb[3L], text)
}

# Position (1 = worst, 10 = best) of a value on the colour ramp for one kind of statistic.
psw_color_index <- function(kind, value) {
  switch(kind,
    cfi = , tli = findInterval(value, c(-Inf, .70, .75, .80, .85, .90, .93, .95, .97, .99)),
    rmsea = 11L - findInterval(value, c(-Inf, .01, .05, .08, .10, .12, .15, .20, .25, .30)),
    srmr = 11L - findInterval(value, c(-Inf, .03, .05, .08, .10, .15, .20, .25, .30, .35)),
    ave = findInterval(value, c(-Inf, .1, .2, .3, .4, .5, .6, .7, .8, .9)),
    omega = findInterval(value, c(-Inf, .50, .55, .60, .65, .70, .80, .85, .90, .95)),
    loading = findInterval(abs(value), c(-Inf, .1, .2, .3, .4, .5, .6, .7, .8, .9)),
    NA_integer_)
}

psw_cell <- function(value, width = 10L, digits = 3L, kind = "", color = FALSE) {
  if (is.character(value)) return(formatC(value, width = width))
  if (is.na(value)) return(formatC("NA", width = width))
  if (kind == "p" && value < 0.001) return(formatC("<.001", width = width))
  txt <- formatC(value, width = width, digits = digits, format = "f")
  idx <- if (color) psw_color_index(kind, value) else NA_integer_
  if (is.na(idx)) txt else psw_paint(txt, idx)
}

psw_short <- function(x, width) {
  ifelse(nchar(x) > width, paste0(substr(x, 1L, width - 2L), ".."), x)
}

# Print a block: first column of row labels, then one column per `cols`.
psw_print_block <- function(title, row_labels, cols, cells, label_width = NULL) {
  if (is.null(label_width)) label_width <- max(8L, min(28L, max(nchar(row_labels)) + 1L))
  cat("\n", title, "\n", sep = "")
  cat(formatC("", width = -label_width), paste(formatC(psw_short(cols, 9L), width = 10L), collapse = ""), "\n", sep = "")
  for (i in seq_along(row_labels))
    cat(formatC(psw_short(row_labels[i], label_width - 1L), width = -label_width),
        paste(cells[[i]], collapse = ""), "\n", sep = "")
}

psw_print_notes <- function(notes) {
  if (!length(notes)) return(invisible())
  cat("\nNotes\n")
  for (n in notes) cat(strwrap(n, width = 0.9 * getOption("width"), initial = "- ", exdent = 2L), sep = "\n")
}

psw_header <- function(title, colored) {
  if (colored) cat("\033[1;36m", title, "\033[0m\n", sep = "") else cat(title, "\n", sep = "")
}

#' Print psyWiz results
#'
#' @param x A psyWiz result: `psw_multi`, `psw_invariance`, `psw_scale`, `psw_esem`,
#'   `psw_multiverse` or `psw_composite`.
#' @param colored Colour-code fit indices and loadings on a red-to-blue ramp
#'   (only when the terminal supports it)? Off by default.
#' @param loadings For `psw_multi`: `"std"`, `"unstd"`, `"both"` or `"none"`.
#' @param sort For `psw_multi`: sort items by mean absolute standardised
#'   loading within each factor?
#' @param culprits For `psw_invariance`: number of worst-fitting constraints
#'   to list for each step that is not supported (0 for none).
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @name print.psyWiz
NULL

#' @rdname print.psyWiz
#' @export
print.psw_multi <- function(x, colored = FALSE, loadings = c("std", "unstd", "both", "none"),
                            sort = FALSE, ...) {
  loadings <- match.arg(loadings)
  col <- psw_use_color(colored)
  est <- x$info$estimators
  psw_header("psyWiz: CFA under several estimators", col)
  cat(sprintf("%d factor(s), %d items, N = %d, items are %s\n", length(x$info$factors),
              length(unlist(x$info$factors)), x$info$n, x$info$type))

  f <- x$fit
  rows <- c("treated as", "N", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.lo", "rmsea.hi", "srmr")
  cells <- lapply(rows, function(r) vapply(seq_along(est), function(k) {
    if (r == "treated as") return(psw_cell(substr(f$treated[k], 1L, 4L)))
    if (r == "N") return(psw_cell(f$n[k], digits = 0L))
    if (r == "df") return(psw_cell(f$df[k], digits = if (isTRUE(f$df[k] == round(f$df[k]))) 0L else 2L))
    psw_cell(f[[r]][k], kind = if (r == "pvalue") "p" else sub("\\..*$", "", r), color = col)
  }, ""))
  psw_print_block("Fit indices", rows, est, cells)

  show <- function(value, title) {
    w <- psw_wide(x$loadings, c("factor", "item"), value, est)
    if (is.null(w)) return(invisible())
    if (sort) {
      avg <- rowMeans(abs(as.matrix(w[est])), na.rm = TRUE)
      w <- w[order(match(w$factor, names(x$info$factors)), -avg), , drop = FALSE]
    }
    lab <- if (length(x$info$factors) > 1L) paste0(psw_short(w$factor, 8L), ": ", w$item) else w$item
    cells <- lapply(seq_len(nrow(w)), function(i) vapply(est, function(e)
      psw_cell(w[[e]][i], kind = if (value == "std") "loading" else "", color = col), ""))
    psw_print_block(title, lab, est, cells)
  }
  if (loadings %in% c("unstd", "both")) show("est", "Unstandardized loadings")
  if (loadings %in% c("std", "both")) show("std", "Standardized loadings")

  w <- psw_wide(x$cors, c("factor1", "factor2"), "r", est)
  if (!is.null(w)) {
    cells <- lapply(seq_len(nrow(w)), function(i) vapply(est, function(e) psw_cell(w[[e]][i]), ""))
    psw_print_block("Factor correlations", paste(w$factor1, "~", w$factor2), est, cells)
  }
  for (v in c("omega", "ave")) {
    w <- psw_wide(x$reliability, "factor", v, est)
    if (is.null(w)) next
    cells <- lapply(seq_len(nrow(w)), function(i) vapply(est, function(e)
      psw_cell(w[[e]][i], kind = v, color = col), ""))
    psw_print_block(if (v == "omega") "Composite reliability (omega)" else "Average variance extracted (AVE)",
                    w$factor, est, cells)
  }
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_invariance <- function(x, colored = FALSE, culprits = 5L, ...) {
  col <- psw_use_color(colored)
  psw_header("psyWiz: measurement invariance across groups", col)
  cat(sprintf("%d factor(s), %d items, items are %s, estimator %s, missing data: %s\n",
              length(x$info$factors), length(unlist(x$info$factors)), x$info$type,
              x$info$estimator, x$info$missing))
  cat("Groups (", x$info$group, "): ",
      paste(sprintf("%s n = %d", x$groups$group, x$groups$n), collapse = ", "), "\n", sep = "")

  f <- x$fit
  rows <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")
  cells <- lapply(seq_len(nrow(f)), function(i) vapply(rows, function(r) {
    if (r == "df") return(psw_cell(f$df[i], digits = if (isTRUE(f$df[i] == round(f$df[i]))) 0L else 2L))
    psw_cell(f[[r]][i], kind = if (r == "pvalue") "p" else r, color = col)
  }, ""))
  psw_print_block("Model fit", f$label, rows, cells)

  cp <- x$comparison
  mark <- function(v, flag, kind = "") {
    txt <- psw_cell(v, kind = kind)
    if (isTRUE(flag)) { if (col) psw_paint(txt, 2L) else sub(" (?=\\S)", "!", txt, perl = TRUE) } else txt
  }
  cells <- lapply(seq_len(nrow(cp)), function(i) c(
    mark(cp$dchisq[i], FALSE), psw_cell(cp$ddf[i], digits = 0L), mark(cp$p[i], cp$flag.chisq[i], "p"),
    mark(cp$dcfi[i], cp$flag.cfi[i]), mark(cp$drmsea[i], cp$flag.rmsea[i]),
    mark(cp$dsrmr[i], cp$flag.srmr[i]), paste0("  ", cp$verdict[i])))
  cat("(each step keeps the constraints of the steps above it)\n")
  psw_print_block("Each step against the previous one",
                  cp$label, c("d.chisq", "d.df", "p", "d.cfi", "d.rmsea", "d.srmr", "  verdict"), cells)
  cut <- x$cutoffs
  cat(sprintf("Verdict uses: %s. Cutoffs: CFI drop > %.3f, RMSEA rise > %.3f, SRMR rise > %.3f (loadings) / %.3f, p < %.2f%s\n",
              paste(x$criteria, collapse = " + "), cut$cfi, cut$rmsea, cut$srmr.loadings, cut$srmr, cut$alpha,
              if (col) "" else "; '!' marks a flagged value"))

  if (culprits > 0L && !is.null(x$culprits)) {
    for (s in cp$step[cp$verdict %in% c("Not supported", "Ambiguous")]) {
      d <- utils::head(x$culprits[x$culprits$step == s, , drop = FALSE], culprits)
      if (!nrow(d)) next
      lab <- paste0(d$item, ifelse(nzchar(d$detail), paste0(" ", d$detail), ""), " [", d$group, "]")
      cells <- lapply(seq_len(nrow(d)), function(i) c(psw_cell(d$x2[i]), psw_cell(d$p[i], kind = "p")))
      psw_print_block(sprintf("Worst-fitting constraints, %s (score tests)", cp$label[cp$step == s]),
                      lab, c("chisq", "p"), cells)
    }
  }
  if (!is.null(x$partial)) {
    cat("\n")
    print(x$partial, colored = colored)
  }
  if (!is.null(x$means)) {
    cat("\n")
    print(x$means, colored = colored)
  }
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_means <- function(x, ...) {
  if (!isTRUE(attr(x, "supported"))) {
    psw_header("psyWiz: latent mean comparisons across groups", FALSE)
    cat("Scalar invariance is not supported across groups. Latent means cannot be compared.\n")
    return(invisible(x))
  }
  ref <- attr(x, "reference")
  mod_type <- attr(x, "model_type")
  psw_header("psyWiz: latent mean comparisons across groups", FALSE)
  cat(sprintf("Reference group: %s (baseline mean = 0.00)%s\n", ref,
              if (identical(mod_type, "partial")) " [partial scalar invariance]" else " [full scalar invariance]"))
  if (identical(mod_type, "partial") && length(attr(x, "freed_items")))
    cat(sprintf("Freed items: %s\n", paste(unique(attr(x, "freed_items")), collapse = ", ")))

  df <- as.data.frame(x)
  fmt_num <- function(val, digits = 3L) {
    ifelse(is.na(val), "      NA", formatC(val, digits = digits, format = "f", width = 8L))
  }
  fmt_p <- function(p) {
    ifelse(is.na(p), "      NA", ifelse(p < 0.001, "   <.001", formatC(p, digits = 3L, format = "f", width = 8L)))
  }
  ci_txt <- sprintf("[%s, %s]", fmt_num(df$ci_lower), fmt_num(df$ci_upper))
  out_df <- data.frame(
    Factor = df$factor,
    Group = df$group,
    Diff = fmt_num(df$mean),
    SE = fmt_num(df$se),
    `95% CI` = ci_txt,
    z = fmt_num(df$z, 2L),
    p = fmt_p(df$p),
    `d (ref)` = fmt_num(df$d_ref),
    `d (pool)` = fmt_num(df$d_pool),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  print(out_df, row.names = FALSE)
  psw_print_notes(attr(x, "notes"))
  invisible(x)
}

#' @export
print.psw_partial <- function(x, colored = FALSE, ...) {
  col <- psw_use_color(colored)
  psw_header("psyWiz: partial measurement invariance", col)

  cp <- x$comparison
  mark <- function(v, flag, kind = "") {
    txt <- psw_cell(v, kind = kind)
    if (isTRUE(flag)) { if (col) psw_paint(txt, 2L) else sub(" (?=\\S)", "!", txt, perl = TRUE) } else txt
  }
  cells <- lapply(seq_len(nrow(cp)), function(i) c(
    mark(cp$dchisq[i], FALSE), psw_cell(cp$ddf[i], digits = 0L), mark(cp$p[i], cp$flag.chisq[i], "p"),
    mark(cp$dcfi[i], cp$flag.cfi[i]), mark(cp$drmsea[i], cp$flag.rmsea[i]),
    mark(cp$dsrmr[i], cp$flag.srmr[i]), paste0("  ", cp$verdict[i]),
    paste0("  ", cp$freed[i])))
  psw_print_block("Partial invariance ladder",
                  cp$label, c("d.chisq", "d.df", "p", "d.cfi", "d.rmsea", "d.srmr", "  verdict", "  freed"), cells)

  if (!is.null(x$released) && nrow(x$released) > 0L) {
    rel <- x$released
    row_labs <- paste0(rel$step, ": ", rel$item, " (", rel$parameter, ")")
    rcells <- lapply(seq_len(nrow(rel)), function(i) c(
      psw_cell(rel$group[i]),
      psw_cell(rel$x2[i]),
      psw_cell(rel$p[i], kind = "p"),
      paste0("  ", rel$reason[i])
    ))
    psw_print_block("Released constraints", row_labs, c("group", "score.x2", "p", "  reason"), rcells)
  }

  if (nzchar(x$summary)) {
    cat("\nSubstantive conclusion:\n  ", x$summary, "\n", sep = "")
  }
  invisible(x)
}

# ---- print methods for the single-scale, ESEM, multiverse and composite results ----

psw_print_df <- function(title, d, digits = 3L) {
  cat("\n", title, "\n", sep = "")
  whole <- vapply(d, is.integer, logical(1))
  d[whole] <- lapply(d[whole], as.character)
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f")))
  d[!num] <- lapply(d[!num], function(v) ifelse(is.na(v), "", as.character(v)))
  print(d, row.names = FALSE, right = FALSE)
}

#' @rdname print.psyWiz
#' @export
print.psw_scale <- function(x, colored = FALSE, ...) {
  psw_header(sprintf("psyWiz: scale check for %s", x$info$name), psw_use_color(colored))
  cat(sprintf("%d items, items are %s, estimator %s, N = %d (%d complete)\n", nrow(x$items), x$info$type,
              x$info$estimator, x$info$n, x$info$n.complete))
  psw_print_df("Reliability", x$reliability)
  psw_print_df("Measurement models", x$fit[c("label", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")])
  psw_print_df("Each model against the previous one",
               x$comparison[c("label", "dchisq", "ddf", "p", "dcfi", "drmsea", "dsrmr", "verdict")])
  cat(paste0("  ", psw_scale_reading(x)), sep = "\n")
  psw_print_df("Items", x$items[c("item", "loading", "r2", "rest.r", "alpha.drop", "H")])
  psw_print_df("Scalability", x$scalability)
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_esem <- function(x, colored = FALSE, ...) {
  psw_header("psyWiz: CFA against ESEM", psw_use_color(colored))
  cat(sprintf("%d factors, %d items, items are %s, estimator %s, N = %d\n", length(x$info$factors),
              length(unlist(x$info$factors)), x$info$type, x$info$estimator, x$info$n))
  psw_print_df("Model fit", x$fit[c("label", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")])
  psw_print_df("CFA minus ESEM", x$comparison[c("dchisq", "ddf", "p", "dcfi", "drmsea", "dsrmr")])
  L <- x$loadings
  fn <- names(x$info$factors)
  items <- unique(L[c("own", "item")])
  wide <- cbind(items, CFA = L$cfa[L$target][match(items$item, L$item[L$target])])
  for (f in fn) wide[[paste0("ESEM:", f)]] <- L$esem[L$factor == f][match(items$item, L$item[L$factor == f])]
  psw_print_df("Standardized loadings", wide)
  psw_print_df("Factor correlations by method", x$correlations)
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_multiverse <- function(x, colored = FALSE, ...) {
  psw_header("psyWiz: robustness of the invariance verdict", psw_use_color(colored))
  s <- x$summary
  psw_print_df("How often each step is supported",
               data.frame(step = s$label, supported = sprintf("%d of %d", s$supported, s$total),
                          share = s$share, reading = s$reading))
  psw_print_df("Verdict under each specification", psw_multiverse_wide(x))
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_composite <- function(x, colored = FALSE, ...) {
  psw_header("psyWiz: is one composite enough?", psw_use_color(colored))
  cat(sprintf("%d indicators, %d outcomes, N = %d\n\n", length(x$info$indicators), length(x$info$outcomes), x$info$n))
  cat(strwrap(x$summary, width = 0.9 * getOption("width"), prefix = "  "), sep = "\n")
  psw_print_df("Rank-one test (Bartlett)", x$test)
  psw_print_df(sprintf("R-squared per outcome (cv: %d folds x %d)", x$info$folds, x$info$repeats), x$r2)
  psw_print_df("Composite weights", x$weights)
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_cfa <- function(x, colored = FALSE, suggestions = 5L, ...) {
  psw_header("psyWiz: confirmatory factor analysis", psw_use_color(colored))
  cat(sprintf("%d factor(s), %d items, items are %s, estimator %s, N = %d\n", length(x$info$factors),
              length(unlist(x$info$factors)), x$info$type, x$info$estimator, x$info$n))
  psw_print_df("Model fit", x$fit[c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.lo", "rmsea.hi", "srmr")])
  psw_print_df("Loadings", x$loadings[c("factor", "item", "est", "se", "p", "std", "r2")])
  if (!is.null(x$cors) && nrow(x$cors)) psw_print_df("Factor correlations", x$cors[c("factor1", "factor2", "r", "p")])
  if (!is.null(x$higher)) psw_print_df(paste("Loadings on", x$higher$name), x$higher$loadings)
  if (!is.null(x$residuals)) psw_print_df("Residual correlations", x$residuals)
  psw_print_df("Reliability", x$reliability)
  if (!is.null(x$correlations)) psw_print_df("Factor correlations by scoring method", x$correlations)
  if (!is.null(x$suggestions) && suggestions > 0L)
    psw_print_df("Largest modification indices", utils::head(x$suggestions, suggestions))
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_bifactor <- function(x, colored = FALSE, ...) {
  psw_header("psyWiz: bifactor CFA", psw_use_color(colored))
  cat(sprintf("%s + %d specific factors, %d items, items are %s, estimator %s, N = %d\n", x$info$general,
              length(x$info$factors), nrow(x$loadings), x$info$type, x$info$estimator, x$info$n))
  psw_print_df("Model fit", x$fit[c("label", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")])
  psw_print_df("Omega, omega hierarchical and explained common variance", x$indices)
  psw_print_df("Standardized loadings", x$loadings)
  psw_print_notes(x$notes)
  invisible(x)
}

#' @rdname print.psyWiz
#' @export
print.psw_order <- function(x, colored = FALSE, patterns = 8L, ...) {
  psw_header("psyWiz: order among items (exploratory)", psw_use_color(colored))
  cat(sprintf("%d items, N = %d, exceptions tolerated %.0f%%, H at least %.2f\n\n", x$info$k, x$info$n,
              100 * x$info$tolerance, x$info$h.min))
  cat(paste0("  ", x$reading), sep = "\n")
  psw_print_df("How well does an order fit?", x$fit)
  psw_print_df("Items", x$items)
  if (!is.null(x$implications)) psw_print_df("Direct implications (from implies to)", x$implications)
  if (!is.null(x$scales)) psw_print_df("Cumulative scales", x$scales)
  if (patterns > 0L) psw_print_df("Most frequent answer patterns", utils::head(x$patterns, patterns))
  psw_print_notes(x$notes)
  invisible(x)
}

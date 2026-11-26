# plot.R -- figures.
#
#   psw_diagram()         path diagram (fixed layout)
#   psw_comb()            comb plot, also for second-order models
#   psw_scalogram_plot()  Guttman scalogram
#   psw_antenna(), psw_facing(), psw_bifactor_path()   bifactor models
#   psw_brush()           loadings of every item on every factor (ESEM)
#   plot() methods at the end
#
# All figures are black, white and grey. Colouring the marks was tried in
# 0.4 (tooth colour by loading size) and dropped: it added nothing that the
# length does not already show.
#
# ---- path diagram ----
#
# General-purpose path-diagram layouts tend to look bad because they have to
# cope with any model. psyWiz only fits simple-structure CFAs (every item on
# one factor, factors correlated), so the diagram uses one fixed layout that
# always works for that kind of model:
#
#   factors in a column on the left, items in a column on the right, one
#   straight arrow per loading, factor correlations as arcs on the far left;
#   one panel per group (or per estimator), all panels on the same layout.
#
# The line width of an arrow is proportional to the standardised loading, so
# panels can be compared at a glance; the number sits on the arrow near the
# item. Nothing is positioned by an optimiser and nothing can overlap.
#
# Text widths are not known to ggplot2 in data units, so the space reserved
# for labels is estimated from the number of characters, assuming the image
# size returned by psw_diagram_size() (about 110 px per x unit, 30 px per row).

psw_dg <- list(char = 0.066, px_x = 110, px_y = 30, max_factor = 22L, max_item = 18L,
               max_cor_factors = 4L, weak = 0.40)

#' Format a single-line summary of model fit (internal)
#'
#' @param fit A data frame or list with fit measures (`chisq`, `df`, `pvalue`, `cfi`, `tli`, `rmsea`, `srmr`), or a character string.
#' @param estimator Name of estimator (e.g. "MLR", "WLSMV").
#' @return A character string, or NULL if no fit measures are available.
#' @keywords internal
#' @noRd
psw_format_fit_line <- function(fit, estimator = NULL) {
  if (is.null(fit)) return(NULL)
  if (is.character(fit) && length(fit) == 1L) return(fit)
  if (is.data.frame(fit)) {
    if (nrow(fit) == 0L) return(NULL)
    row <- fit[1L, , drop = FALSE]
  } else if (is.list(fit)) {
    row <- as.data.frame(fit)
  } else {
    return(NULL)
  }

  parts <- character()
  if (!is.null(estimator) && is.character(estimator) && nzchar(estimator)) {
    parts <- c(parts, paste0("Estimator: ", estimator))
  }

  fmt_num <- function(x, digits = 3, lead_zero = TRUE) {
    if (is.null(x) || is.na(x) || !is.numeric(x)) return(NA_character_)
    s <- sprintf(paste0("%.", digits, "f"), x)
    if (!lead_zero) sub("^0\\.", ".", sub("^-0\\.", "-.", s)) else s
  }

  if ("chisq" %in% names(row) && !is.na(row$chisq)) {
    ch_str <- fmt_num(row$chisq, 2)
    df_str <- if ("df" %in% names(row) && !is.na(row$df)) paste0("(", round(row$df), ")") else ""
    p_str <- if ("pvalue" %in% names(row) && !is.na(row$pvalue)) {
      if (row$pvalue < .001) ", p < .001" else paste0(", p = ", fmt_num(row$pvalue, 3, FALSE))
    } else ""
    parts <- c(parts, paste0("\u03c7\u00b2", df_str, " = ", ch_str, p_str))
  }

  if ("cfi" %in% names(row) && !is.na(row$cfi)) {
    parts <- c(parts, paste0("CFI = ", fmt_num(row$cfi, 3, FALSE)))
  }
  if ("tli" %in% names(row) && !is.na(row$tli)) {
    parts <- c(parts, paste0("TLI = ", fmt_num(row$tli, 3, FALSE)))
  }
  if ("rmsea" %in% names(row) && !is.na(row$rmsea)) {
    ci_str <- ""
    if ("rmsea.lo" %in% names(row) && "rmsea.hi" %in% names(row) &&
        !is.na(row$rmsea.lo) && !is.na(row$rmsea.hi)) {
      ci_str <- paste0(" [", fmt_num(row$rmsea.lo, 3, FALSE), ", ", fmt_num(row$rmsea.hi, 3, FALSE), "]")
    }
    parts <- c(parts, paste0("RMSEA = ", fmt_num(row$rmsea, 3, FALSE), ci_str))
  }
  if ("srmr" %in% names(row) && !is.na(row$srmr)) {
    parts <- c(parts, paste0("SRMR = ", fmt_num(row$srmr, 3, FALSE)))
  }

  if (!length(parts)) return(NULL)
  paste(parts, collapse = " | ")
}

psw_wrap_fit_line <- function(text, width_px = NULL) {
  if (is.null(text) || !nzchar(text)) return(text)
  max_chars <- if (!is.null(width_px) && is.finite(width_px)) max(46L, floor((width_px - 32) / 6.1)) else 72L
  if (nchar(text) <= max_chars) return(text)
  parts <- strsplit(text, " \\| ")[[1L]]
  lines <- character()
  cur <- ""
  for (pt in parts) {
    cand <- if (!nzchar(cur)) pt else paste0(cur, " | ", pt)
    if (nzchar(cur) && nchar(cand) > max_chars) {
      lines <- c(lines, cur)
      cur <- pt
    } else {
      cur <- cand
    }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  paste(lines, collapse = "\n")
}

psw_attach_fit <- function(p, fit_text = NULL, fit_pos = c("none", "below", "above"), width_px = NULL) {
  if (is.null(fit_pos)) fit_pos <- "none"
  fit_pos <- match.arg(fit_pos, c("none", "below", "above"))
  if (fit_pos == "none" || is.null(fit_text) || !nzchar(fit_text)) return(p)

  sz <- attr(p, "psw_size")
  w <- if (!is.null(width_px)) width_px else if (!is.null(sz) && "width" %in% names(sz)) sz[["width"]] else NULL
  wrapped <- psw_wrap_fit_line(fit_text, w)
  n_lines <- length(strsplit(wrapped, "\n", fixed = TRUE)[[1L]])

  if (fit_pos == "above") {
    cur <- p$labels$subtitle
    new_sub <- if (is.null(cur) || !nzchar(cur)) wrapped else paste0(cur, "\n", wrapped)
    p <- p + ggplot2::labs(subtitle = new_sub) +
      ggplot2::theme(plot.title.position = "plot",
                     plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0,
                                                           lineheight = 1.25, margin = ggplot2::margin(2, 0, 8, 0)))
  } else if (fit_pos == "below") {
    cur <- p$labels$caption
    new_cap <- if (is.null(cur) || !nzchar(cur)) wrapped else paste0(wrapped, "\n", cur)
    p <- p + ggplot2::labs(caption = new_cap) +
      ggplot2::theme(plot.caption.position = "plot",
                     plot.caption = ggplot2::element_text(size = 8.8, colour = "grey35", hjust = 0,
                                                          lineheight = 1.25, margin = ggplot2::margin(8, 0, 0, 0)))
  }
  if (!is.null(sz)) {
    sz["height"] <- sz["height"] + 18L * n_lines + 8L
    attr(p, "psw_size") <- sz
  }
  p
}

psw_dg_layout <- function(loadings, cors, panel_col) {
  items <- unique(loadings[c("factor", "item")])
  fnames <- unique(items$factor)
  fi <- match(items$factor, fnames)
  # rows from the top, with a gap between the item blocks of different factors
  items$y <- -(seq_len(nrow(items)) + (fi - 1L) * 0.7)
  fy <- tapply(items$y, factor(items$factor, levels = fnames), mean)
  flab <- psw_short(fnames, psw_dg$max_factor)
  ilab <- psw_short(items$item, psw_dg$max_item)
  wf <- max(nchar(flab)) * psw_dg$char + 0.08
  wi <- max(nchar(ilab)) * psw_dg$char + 0.08
  show_cors <- !is.null(cors) && nrow(cors) > 0L && length(fnames) <= psw_dg$max_cor_factors
  span <- if (show_cors) max(abs(outer(fy, fy, "-"))) else 0
  # room for the arcs grows with the longest arc (in rows)
  arc <- if (show_cors) 0.22 + 0.035 * span else 0
  list(items = items, fnames = fnames, fy = fy, flab = flab, ilab = ilab, wf = wf, wi = wi,
       show_cors = show_cors, arc = arc,
       xlim = c(-(wf + arc + if (show_cors) 0.52 else 0.02), 1 + wi + 0.04),
       ylim = c(min(items$y) - 0.7, -0.3),
       panels = unique(loadings[[panel_col]]))
}

# Image size in pixels that matches the layout assumptions above.
psw_diagram_size <- function(loadings, cors = NULL, panel_col = "group") {
  lay <- psw_dg_layout(loadings, cors, panel_col)
  np <- length(lay$panels)
  ncol <- if (np <= 4L) np else 3L
  nrow <- ceiling(np / ncol)
  c(width = round(ncol * diff(lay$xlim) * psw_dg$px_x + 30),
    height = round(nrow * (diff(lay$ylim) * psw_dg$px_y + 30) + 50))
}

psw_fmt_loading <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))

#' Draw a simple-structure CFA as a fixed-layout diagram (internal)
#'
#' @param loadings Data frame with columns `factor`, `item`, `std` and the
#'   panel column.
#' @param cors Data frame with `factor1`, `factor2`, `r` and the panel column,
#'   or `NULL`.
#' @param panel_col Name of the column that defines the panels.
#' @param numbers Print the standardized loadings on the arrows?
#' @param width Make the width of an arrow proportional to its loading?
#' @return A ggplot object.
#' @keywords internal
#' @noRd
psw_diagram <- function(loadings, cors = NULL, panel_col = "group", caption = NULL, numbers = TRUE,
                        width = TRUE, fit_text = NULL, fit_pos = c("none", "below", "above")) {
  lay <- psw_dg_layout(loadings, cors, panel_col)
  it <- lay$items
  panel <- function(d) factor(d[[panel_col]], levels = lay$panels)
  ink <- "grey25"

  m <- match(paste(loadings$factor, loadings$item), paste(it$factor, it$item))
  arrows <- data.frame(panel = panel(loadings), x = 0.03, xend = 0.97,
                       y = unname(lay$fy[loadings$factor]), yend = it$y[m],
                       std = loadings$std, stringsAsFactors = FALSE)
  arrows <- arrows[is.finite(arrows$std), , drop = FALSE]
  arrows$lw <- if (width) 0.25 + 1.5 * pmin(abs(arrows$std), 1) else 0.5
  arrows$col <- ink
  # the number sits on the arrow, near the item, where the fan has spread out
  t <- 0.68
  arrows$lx <- arrows$x + t * (arrows$xend - arrows$x)
  arrows$ly <- arrows$y + t * (arrows$yend - arrows$y)
  arrows$lab <- psw_fmt_loading(arrows$std)

  grid <- expand.grid(i = seq_len(nrow(it)), panel = factor(lay$panels, levels = lay$panels))
  inodes <- data.frame(panel = grid$panel, x = 1, y = it$y[grid$i], lab = lay$ilab[grid$i])
  gridf <- expand.grid(i = seq_along(lay$fnames), panel = factor(lay$panels, levels = lay$panels))
  fnodes <- data.frame(panel = gridf$panel, x = 0, y = unname(lay$fy[gridf$i]), lab = lay$flab[gridf$i])

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                                  linewidth = .data$lw, colour = .data$col),
      lineend = "round", arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches"), type = "closed"))
  if (numbers) p <- p + ggplot2::geom_label(
      data = arrows, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$lab, colour = .data$col),
      size = 3.1, fill = "white", linewidth = 0, label.padding = ggplot2::unit(0.08, "lines"))
  p <- p + ggplot2::geom_label(
    data = inodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
    hjust = 0, size = 3.5, colour = ink, fill = "white", linewidth = 0.3,
    label.r = ggplot2::unit(0, "lines"), label.padding = ggplot2::unit(0.3, "lines"))
  p <- p + ggplot2::geom_label(
    data = fnodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
    hjust = 1, size = 3.7, fontface = "bold", colour = ink, fill = "grey93", linewidth = 0.3,
    label.r = ggplot2::unit(0.6, "lines"), label.padding = ggplot2::unit(0.45, "lines"))

  if (lay$show_cors) {
    x0 <- -(lay$wf + 0.09)
    arcs <- data.frame(panel = panel(cors), x = x0,
                       y = unname(lay$fy[cors$factor1]), yend = unname(lay$fy[cors$factor2]),
                       lab = psw_fmt_loading(cors$r), stringsAsFactors = FALSE)
    arcs$lx <- x0 - lay$arc - 0.03 * abs(arcs$y - arcs$yend)
    arcs$ly <- (arcs$y + arcs$yend) / 2
    p <- p +
      ggplot2::geom_curve(
        data = arcs, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$x, yend = .data$yend),
        curvature = 0.55, linewidth = 0.3, colour = "grey45", ncp = 12,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches"), ends = "both", type = "closed")) +
      ggplot2::geom_text(data = arcs, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$lab),
                         size = 3.1, colour = "grey35", hjust = 1)
  }

  ncol <- if (length(lay$panels) <= 4L) length(lay$panels) else 3L
  p <- p +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_wrap(~panel, ncol = ncol) +
    ggplot2::coord_cartesian(xlim = lay$xlim, ylim = lay$ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(caption = caption) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 12,
                                                      margin = ggplot2::margin(4, 0, 6, 0)),
                   plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0,
                                                        margin = ggplot2::margin(8, 0, 0, 0)),
                   panel.spacing = ggplot2::unit(14, "pt"),
                   plot.margin = ggplot2::margin(8, 10, 8, 10),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  psw_attach_fit(p, fit_text, fit_pos)
}

# ---- comb plot ---------------------------------------------------------------
#
# One comb per factor. The spine is the factor, with its name written along
# it. Every item is a tooth whose dark part has length loading^2, the share of
# the item's variance it has in common with the factor, on a pale tooth of
# length 1 (so the pale remainder is the item's unique variance). Every
# length is on the same 0-1 scale, which is what makes combs from different
# groups, estimators or scales directly comparable.
#
#   factor correlations    arcs joining the spines, thicker = stronger
#                          (up to four factors; beyond that they would tangle)
#   residual correlations  brackets joining two tooth tips (only those in the model)
#   thick = TRUE           tooth thickness proportional to the item's observed variance
#   higher                 second-order model: each factor gets a top tooth whose dark
#                          part is the share of the factor explained by the higher
#                          factor, and one more spine on the left joins those teeth
#
# History: up to 0.3 the factor name sat on a top tooth of length 1 and the
# correlations were thin teeth under it. The variants were compared with
# tools/plot_proposals.R on 2026-11-16 and the author chose this one.

psw_cb <- list(px_x = 210, px_y = 22, char = 7.2)

psw_comb_layout <- function(loadings, cors, panel_col, resid = NULL, higher = NULL) {
  items <- unique(loadings[c("factor", "item")])
  fnames <- unique(items$factor)
  nf <- length(fnames)
  show_arcs <- is.null(higher) && !is.null(cors) && nrow(cors) > 0L && nf >= 2L && nf <= psw_dg$max_cor_factors
  y <- 0
  top <- bottom <- stats::setNames(numeric(nf), fnames)
  items$y <- NA_real_
  for (f in fnames) {
    idx <- which(items$factor == f)
    top[f] <- y
    items$y[idx] <- y - seq_along(idx) + if (is.null(higher)) 1 else -0.35
    bottom[f] <- min(items$y[idx])
    y <- bottom[f] - 2.4
  }
  ilab <- psw_short(items$item, 22L)
  wl <- (max(nchar(ilab)) * psw_cb$char + 14) / psw_cb$px_x
  # the name runs along the spine, so it can be as long as the spine is
  room <- pmax(5, floor((top - bottom + 1) * psw_cb$px_y / 7.8))
  flab <- vapply(seq_len(nf), function(i) psw_short(fnames[i], room[[i]]), "")
  name_x <- -(wl + 0.07)
  arc_x <- name_x - 0.08
  depth <- function(span) 0.15 + 0.12 * (span - 1)
  left <- if (show_arcs) arc_x - depth(nf - 1L) - 0.16 else if (!is.null(higher)) arc_x - 0.3 else name_x - 0.08
  right <- if (!is.null(resid) && nrow(resid) > 0L) 0.5 else 0.17
  list(items = items, ilab = ilab, fnames = fnames, flab = flab, top = top, bottom = bottom,
       show_arcs = show_arcs, name_x = name_x, arc_x = arc_x, depth = depth,
       xlim = c(left, 1 + right), ylim = c(min(items$y) - 0.8, 0.9),
       panels = unique(loadings[[panel_col]]))
}

psw_comb_size <- function(loadings, cors = NULL, panel_col = "group", resid = NULL, higher = NULL) {
  lay <- psw_comb_layout(loadings, cors, panel_col, resid, higher)
  np <- length(lay$panels)
  ncol <- if (np <= 3L) np else 3L
  nrow <- ceiling(np / ncol)
  c(width = max(520, round(ncol * diff(lay$xlim) * psw_cb$px_x + 30)),
    height = round(nrow * (diff(lay$ylim) * psw_cb$px_y + 28) + 55))
}

#' Comb plot of a simple-structure factor model (internal)
#' @inheritParams psw_diagram
#' @param thick Tooth thickness proportional to the observed variance of the
#'   item (column `var` of `loadings`)?
#' @param resid Data frame of residual correlations in the model (`item1`,
#'   `item2`, `r`, and the panel column), or `NULL`.
#' @param higher Second-order model: `list(name, loadings = data.frame(factor, std))`.
#' @keywords internal
#' @noRd
psw_comb <- function(loadings, cors = NULL, panel_col = "group", caption = NULL, thick = FALSE,
                     resid = NULL, higher = NULL, fit_text = NULL, fit_pos = c("none", "below", "above")) {
  lay <- psw_comb_layout(loadings, cors, panel_col, resid, higher)
  it <- lay$items
  lv <- lay$panels
  ink <- "grey20"
  panel <- function(d) factor(d[[panel_col]], levels = lv)
  pg <- function(n) factor(rep(lv, each = n), levels = lv)
  nf <- length(lay$fnames)

  m <- match(paste(loadings$factor, loadings$item), paste(it$factor, it$item))
  teeth <- data.frame(panel = panel(loadings), y = it$y[m], std = loadings$std, stringsAsFactors = FALSE)
  v <- if (thick && !is.null(loadings$var) && any(is.finite(loadings$var))) loadings$var else NULL
  thick <- !is.null(v)
  teeth$lw <- if (thick) 0.9 + 3.6 * v / max(v, na.rm = TRUE) else 2.6
  teeth$lw[!is.finite(teeth$lw)] <- 2.6
  teeth <- teeth[is.finite(teeth$std), , drop = FALSE]
  teeth$len <- pmin(teeth$std^2, 1)
  teeth$lab <- psw_fmt_loading(teeth$std)

  inames <- data.frame(panel = pg(nrow(it)), y = rep(it$y, length(lv)), lab = rep(lay$ilab, length(lv)))
  spine <- data.frame(panel = pg(nf), y = rep(unname(lay$top), length(lv)), yend = rep(unname(lay$bottom), length(lv)),
                      lab = rep(lay$flab, length(lv)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = teeth, ggplot2::aes(x = 0, xend = 1, y = .data$y, yend = .data$y,
                                                     linewidth = .data$lw), colour = "grey88") +
    ggplot2::geom_segment(data = teeth, ggplot2::aes(x = 0, xend = .data$len, y = .data$y, yend = .data$y,
                                                     linewidth = .data$lw), colour = ink) +
    ggplot2::geom_text(data = teeth, ggplot2::aes(x = .data$len + 0.02, y = .data$y, label = .data$lab),
                       hjust = 0, size = 3.2, colour = ink) +
    ggplot2::geom_text(data = inames, ggplot2::aes(x = -0.035, y = .data$y, label = .data$lab),
                       hjust = 1, size = 3.4, colour = ink) +
    ggplot2::geom_text(data = spine, ggplot2::aes(x = lay$name_x, y = (.data$y + .data$yend) / 2, label = .data$lab),
                       angle = 90, size = 3.9, fontface = "bold", colour = ink)

  if (lay$show_arcs) {
    mid <- (lay$top + lay$bottom) / 2
    t <- seq(0, pi, length.out = 50L)
    arcs <- do.call(rbind, lapply(seq_len(nrow(cors)), function(i) {
      if (!is.finite(cors$r[i])) return(NULL)
      y1 <- mid[[cors$factor1[i]]]; y2 <- mid[[cors$factor2[i]]]
      d <- lay$depth(abs(match(cors$factor1[i], lay$fnames) - match(cors$factor2[i], lay$fnames)))
      data.frame(panel = panel(cors)[i], id = i, x = lay$arc_x - d * sin(t), y = y1 + (y2 - y1) * (1 - cos(t)) / 2,
                 lw = 0.4 + 3.6 * min(abs(cors$r[i]), 1), lx = lay$arc_x - d, ly = (y1 + y2) / 2,
                 lab = psw_fmt_loading(cors$r[i]), stringsAsFactors = FALSE)
    }))
    if (!is.null(arcs)) {
      labs <- arcs[!duplicated(arcs$id), , drop = FALSE]
      p <- p +
        ggplot2::geom_path(data = arcs, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id, linewidth = .data$lw),
                           colour = "grey40", lineend = "round") +
        ggplot2::geom_label(data = labs, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$lab),
                            size = 3, colour = ink, fill = "white", label.padding = ggplot2::unit(0.12, "lines"))
    }
  }

  if (!is.null(higher)) {
    h <- higher$loadings
    tt <- data.frame(panel = pg(nf), y = rep(unname(lay$top), length(lv)),
                     std = rep(h$std[match(lay$fnames, h$factor)], length(lv)))
    tt$len <- pmin(tt$std^2, 1)
    hx <- lay$arc_x - 0.12
    p <- p +
      ggplot2::geom_segment(data = tt, ggplot2::aes(x = 0, xend = 1, y = .data$y, yend = .data$y), colour = "grey70", linewidth = 3.4) +
      ggplot2::geom_segment(data = tt, ggplot2::aes(x = 0, xend = .data$len, y = .data$y, yend = .data$y), colour = ink, linewidth = 3.4) +
      ggplot2::geom_text(data = tt, ggplot2::aes(x = .data$len + 0.02, y = .data$y + 0.45, label = psw_fmt_loading(.data$std)),
                         hjust = 0, size = 3.1, colour = ink) +
      ggplot2::geom_segment(data = tt, ggplot2::aes(x = hx, xend = 0, y = .data$y, yend = .data$y), colour = ink,
                            linewidth = 0.5, linetype = "dotted") +
      ggplot2::annotate("segment", x = hx, xend = hx, y = max(lay$top), yend = min(lay$top), colour = ink,
                        linewidth = 3.4, lineend = "square") +
      ggplot2::annotate("text", x = hx - 0.08, y = mean(range(lay$top)), label = psw_short(higher$name, 24L),
                        angle = 90, size = 3.9, fontface = "bold", colour = ink)
  }

  if (!is.null(resid) && nrow(resid) > 0L) {
    br <- data.frame(panel = panel(resid), y1 = it$y[match(resid$item1, it$item)], y2 = it$y[match(resid$item2, it$item)],
                     lab = paste("r =", psw_fmt_loading(resid$r)), stringsAsFactors = FALSE)
    br <- br[is.finite(br$y1) & is.finite(br$y2), , drop = FALSE]
    # brackets that overlap are pushed outwards so each can be followed
    br <- br[order(abs(br$y1 - br$y2)), , drop = FALSE]
    br$x <- 1.16 + 0.035 * (seq_len(nrow(br)) - 1L) %% 4L
    path <- do.call(rbind, lapply(seq_len(nrow(br)), function(i)
      data.frame(panel = br$panel[i], id = i, x = c(1.12, br$x[i], br$x[i], 1.12), y = c(br$y1[i], br$y1[i], br$y2[i], br$y2[i]))))
    p <- p +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id), colour = ink, linewidth = 0.6) +
      ggplot2::geom_text(data = br, ggplot2::aes(x = .data$x + 0.02, y = (.data$y1 + .data$y2) / 2, label = .data$lab),
                         hjust = 0, size = 3, colour = ink)
  }

  ncol <- if (length(lv) <= 3L) length(lv) else 3L
  p <- p +
    ggplot2::geom_segment(data = spine, ggplot2::aes(x = 0, xend = 0, y = .data$y, yend = .data$yend),
                          colour = ink, linewidth = 3.4, lineend = "square") +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::facet_wrap(~panel, ncol = ncol) +
    ggplot2::coord_cartesian(xlim = lay$xlim, ylim = lay$ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, caption = caption) +
    psw_theme_bars()
  psw_attach_fit(p, fit_text, fit_pos)
}

# Shared look of the figures made of horizontal bars.
psw_theme_bars <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
                   axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 8.5, colour = "grey40"),
                   axis.title.x = ggplot2::element_text(size = 9, colour = "grey40"),
                   strip.text = ggplot2::element_text(face = "bold", size = 12),
                   plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0,
                                                        margin = ggplot2::margin(8, 0, 0, 0)),
                   panel.spacing = ggplot2::unit(16, "pt"), plot.margin = ggplot2::margin(8, 12, 8, 12),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# One entry point for the two figures of a simple-structure model. `opts` may
# hold thick, resid, higher (comb) and numbers, width (path), fit_text, fit_pos.
psw_figure <- function(loadings, cors, panel_col, type = c("comb", "path"), caption = NULL, opts = list()) {
  type <- match.arg(type)
  fit_text <- opts$fit_text
  fit_pos <- if (!is.null(opts$fit_pos)) opts$fit_pos else "none"
  if (type == "comb") {
    p <- psw_comb(loadings, cors, panel_col, caption, thick = isTRUE(opts$thick), resid = opts$resid, higher = opts$higher)
    sz <- psw_comb_size(loadings, cors, panel_col, opts$resid, opts$higher)
  } else {
    p <- psw_diagram(loadings, cors, panel_col, caption, numbers = !isFALSE(opts$numbers), width = !isFALSE(opts$width))
    sz <- psw_diagram_size(loadings, cors, panel_col)
  }
  attr(p, "psw_size") <- sz
  psw_attach_fit(p, fit_text, fit_pos)
}

# ---- scalogram ---------------------------------------------------------------

#' Scalogram (internal): respondents sorted by total score against items
#'
#' Each row is a group of respondents with similar total scores (low at the
#' bottom), each column an item. Dark = endorsed (binary items) or a high
#' response (polytomous items). In a perfect Guttman scale the dark area is a
#' staircase with nothing dark outside it; the line shows where that staircase
#' would be for binary items.
#' @param sc Output of `psw_scalogram()`.
#' @keywords internal
#' @noRd
psw_scalogram_plot <- function(sc) {
  k <- ncol(sc$z)
  nb <- nrow(sc$z)
  d <- data.frame(x = rep(seq_len(k), each = nb), y = rep(seq_len(nb), k), z = as.numeric(sc$z))
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_raster(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$z)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "grey15"), limits = c(0, 1), breaks = c(0, 0.5, 1),
                                  labels = if (sc$binary) c("0%", "50%", "100%") else c("lowest", "", "highest"),
                                  name = if (sc$binary) "Endorsing" else "Response")
  if (sc$binary) {
    # a group with mean total score s is expected to endorse the s easiest items;
    # white under black so the line shows on dark and on light cells
    st <- data.frame(y = c(seq_len(nb) - 0.5, nb + 0.5), x = c(rowSums(sc$z), sum(sc$z[nb, ])) + 0.5)
    p <- p + ggplot2::geom_step(data = st, ggplot2::aes(x = .data$x, y = .data$y), direction = "vh", colour = "white", linewidth = 1.8) +
      ggplot2::geom_step(data = st, ggplot2::aes(x = .data$x, y = .data$y), direction = "vh", colour = "black", linewidth = 0.6)
  }
  p <- p +
    ggplot2::scale_x_continuous(breaks = seq_len(k), labels = psw_short(sc$items, 14L), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(1, nb), labels = c("lowest", "highest"), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Total score") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = if (k > 6L) 45 else 0, hjust = if (k > 6L) 1 else 0.5),
                   plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  attr(p, "psw_size") <- c(width = max(460, 56 * k + 190), height = 400)
  p
}

# ---- bifactor figures --------------------------------------------------------
# All three take the `loadings` data frame of bifactorWiz(): item, specific
# (name of the item's specific factor), general, own (standardized loadings).

psw_bf_rows <- function(L) {
  sp <- unique(L$specific)
  L$y <- -(seq_len(nrow(L)) + (match(L$specific, sp) - 1L) * 1.2)
  L$g2 <- pmin(L$general^2, 1)
  L$s2 <- pmin(L$own^2, 1)
  blocks <- do.call(rbind, lapply(sp, function(f)
    data.frame(specific = f, y = max(L$y[L$specific == f]) + 0.4, yend = min(L$y[L$specific == f]) - 0.4,
               stringsAsFactors = FALSE)))
  list(L = L, blocks = blocks, wl = (max(nchar(psw_short(L$item, 22L))) * psw_cb$char + 14) / psw_cb$px_x)
}

#' Antenna plot (internal): general factor to the left, specific to the right
#'
#' One row per item. From the mast in the middle, the bar to the left is the
#' share of the item's variance due to the general factor and the bar to the
#' right the share due to its specific factor (squared standardized loadings).
#' An item that leans left is mostly a measure of the general factor.
#' @keywords internal
#' @noRd
psw_antenna <- function(L, general = "General", fit_text = NULL, fit_pos = c("none", "below", "above")) {
  r <- psw_bf_rows(L); L <- r$L
  ink <- "grey20"
  ytop <- max(L$y)
  p <- ggplot2::ggplot(L) +
    ggplot2::geom_segment(ggplot2::aes(x = -1, xend = 1, y = .data$y, yend = .data$y), colour = "grey93", linewidth = 2.6) +
    ggplot2::geom_segment(ggplot2::aes(x = -.data$g2, xend = 0, y = .data$y, yend = .data$y), colour = ink, linewidth = 2.6) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$s2, y = .data$y, yend = .data$y), colour = "grey50", linewidth = 2.6) +
    ggplot2::geom_text(ggplot2::aes(x = -.data$g2 - 0.02, y = .data$y, label = psw_fmt_loading(.data$general)), hjust = 1, size = 3.1, colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = .data$s2 + 0.02, y = .data$y, label = psw_fmt_loading(.data$own)), hjust = 0, size = 3.1, colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = -1.04, y = .data$y, label = psw_short(.data$item, 22L)), hjust = 1, size = 3.4, colour = ink) +
    ggplot2::annotate("segment", x = 0, xend = 0, y = ytop + 0.9, yend = min(L$y) - 0.5, colour = ink, linewidth = 1.6) +
    ggplot2::geom_segment(data = r$blocks, ggplot2::aes(x = 1.04, xend = 1.04, y = .data$y, yend = .data$yend), colour = "grey50", linewidth = 1.2) +
    ggplot2::geom_text(data = r$blocks, ggplot2::aes(x = 1.08, y = (.data$y + .data$yend) / 2, label = psw_short(.data$specific, 16L)),
                       hjust = 0, fontface = "bold", size = 3.6, colour = ink) +
    ggplot2::annotate("text", x = c(-0.5, 0.5), y = ytop + 1.15, label = c(psw_short(general, 20L), "Specific factor"),
                      fontface = "bold", size = 3.7, colour = ink) +
    ggplot2::scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("1", ".50", "0", ".50", "1")) +
    ggplot2::coord_cartesian(xlim = c(-1.06 - r$wl, 1.5), ylim = c(min(L$y) - 0.8, ytop + 1.7), expand = FALSE, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    psw_theme_bars()
  attr(p, "psw_size") <- c(width = round((2.56 + r$wl) * 190 + 30), height = round((ytop - min(L$y) + 2.5) * 24 + 55))
  psw_attach_fit(p, fit_text, fit_pos)
}

#' Facing combs (internal): the bifactor model as two combs sharing their teeth
#'
#' Every row is one item's whole variance, from 0 to 1. The general factor's
#' comb grows from the left spine, the specific factors' combs from the right
#' spines, and what is left in the middle is the item's unique variance, where
#' its name is written. Two combs that nearly touch mean a well explained item.
#' @keywords internal
#' @noRd
psw_facing <- function(L, general = "General", fit_text = NULL, fit_pos = c("none", "below", "above")) {
  r <- psw_bf_rows(L); L <- r$L
  ink <- "grey20"
  L$mid <- (L$g2 + 1 - L$s2) / 2
  p <- ggplot2::ggplot(L) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = .data$y, yend = .data$y), colour = "grey90", linewidth = 4.6) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$g2, y = .data$y, yend = .data$y), colour = ink, linewidth = 4.6) +
    ggplot2::geom_segment(ggplot2::aes(x = 1 - .data$s2, xend = 1, y = .data$y, yend = .data$y), colour = "grey50", linewidth = 4.6) +
    ggplot2::geom_label(ggplot2::aes(x = .data$mid, y = .data$y, label = psw_short(.data$item, 14L)), size = 3, colour = ink,
                        fill = "white", label.padding = ggplot2::unit(0.1, "lines")) +
    ggplot2::geom_text(ggplot2::aes(x = -0.03, y = .data$y, label = psw_fmt_loading(.data$general)), hjust = 1, size = 3.1, colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = 1.03, y = .data$y, label = psw_fmt_loading(.data$own)), hjust = 0, size = 3.1, colour = ink) +
    ggplot2::annotate("segment", x = 0, xend = 0, y = max(L$y) + 0.4, yend = min(L$y) - 0.4, colour = ink, linewidth = 3.4, lineend = "square") +
    ggplot2::annotate("text", x = -0.19, y = mean(range(L$y)), label = psw_short(general, 24L), angle = 90, fontface = "bold", size = 3.9, colour = ink) +
    ggplot2::geom_segment(data = r$blocks, ggplot2::aes(x = 1, xend = 1, y = .data$y, yend = .data$yend), colour = "grey50",
                          linewidth = 3.4, lineend = "square") +
    ggplot2::geom_text(data = r$blocks, ggplot2::aes(x = 1.19, y = (.data$y + .data$yend) / 2,
                                                     label = psw_short(.data$specific, pmax(5, floor((.data$y - .data$yend) * 24 / 7.8)))),
                       angle = 270, fontface = "bold", size = 3.9, colour = ink) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::coord_cartesian(xlim = c(-0.27, 1.27), ylim = c(min(L$y) - 0.8, max(L$y) + 0.8), expand = FALSE, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    psw_theme_bars()
  attr(p, "psw_size") <- c(width = 620, height = round((max(L$y) - min(L$y) + 1.6) * 26 + 60))
  psw_attach_fit(p, fit_text, fit_pos)
}

# Path diagram of a bifactor model: general factor on the left, items in the
# middle, specific factors on the right.
psw_bifactor_path <- function(L, general = "General", numbers = TRUE, width = TRUE,
                              fit_text = NULL, fit_pos = c("none", "below", "above")) {
  r <- psw_bf_rows(L); L <- r$L
  ink <- "grey25"
  wi <- max(nchar(psw_short(L$item, 18L))) * psw_dg$char + 0.1
  sy <- tapply(L$y, factor(L$specific, levels = unique(L$specific)), mean)
  lw <- function(v) if (width) 0.25 + 1.5 * pmin(abs(v), 1) else 0.5
  A <- rbind(data.frame(x = 0.03, y = mean(range(L$y)), xend = 1 - wi / 2 - 0.03, yend = L$y, std = L$general),
             data.frame(x = 2 - 0.03, y = unname(sy[L$specific]), xend = 1 + wi / 2 + 0.03, yend = L$y, std = L$own))
  A$lw <- lw(A$std)
  A$lx <- A$x + 0.7 * (A$xend - A$x); A$ly <- A$y + 0.7 * (A$yend - A$y)
  fn <- data.frame(x = c(0, rep(2, length(sy))), y = c(mean(range(L$y)), unname(sy)),
                   lab = psw_short(c(general, names(sy)), psw_dg$max_factor), hjust = c(1, rep(0, length(sy))))
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = A, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linewidth = .data$lw),
                          colour = ink, lineend = "round", arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches"), type = "closed"))
  if (numbers) p <- p + ggplot2::geom_label(data = A, ggplot2::aes(x = .data$lx, y = .data$ly, label = psw_fmt_loading(.data$std)),
                                            size = 3, colour = ink, fill = "white", linewidth = 0, label.padding = ggplot2::unit(0.08, "lines"))
  wf <- max(nchar(fn$lab)) * psw_dg$char + 0.12
  p <- p +
    ggplot2::geom_label(data = L, ggplot2::aes(x = 1, y = .data$y, label = psw_short(.data$item, 18L)), size = 3.5, colour = ink,
                        fill = "white", linewidth = 0.3, label.r = ggplot2::unit(0, "lines"), label.padding = ggplot2::unit(0.3, "lines")) +
    ggplot2::geom_label(data = fn, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab, hjust = .data$hjust), size = 3.7,
                        fontface = "bold", colour = ink, fill = "grey93", linewidth = 0.3, label.r = ggplot2::unit(0.6, "lines"),
                        label.padding = ggplot2::unit(0.45, "lines")) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::coord_cartesian(xlim = c(-wf, 2 + wf), ylim = c(min(L$y) - 0.7, max(L$y) + 0.7), expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0, margin = ggplot2::margin(8, 0, 0, 0)),
                   plot.margin = ggplot2::margin(8, 10, 8, 10), plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  attr(p, "psw_size") <- c(width = round((2 + 2 * wf) * 150 + 30), height = round((max(L$y) - min(L$y) + 1.4) * psw_dg$px_y + 40))
  psw_attach_fit(p, fit_text, fit_pos)
}

# ---- brush plot --------------------------------------------------------------

#' Brush plot (internal): every item on every factor
#'
#' For models where items load on several factors (EFA, ESEM). One vertical
#' line per factor (the handle) and one row per item; the bar to the right of
#' a factor's line is the squared standardized loading of the item on that
#' factor (a bristle): dark when it is the item's own factor, grey for a
#' cross-loading. A last column shows the item's unique variance as a pale
#' bar. With clean data every factor looks like a toothbrush: a block of
#' long, even bristles on its own items and a bare handle elsewhere. Stray
#' bristles are cross-loadings; short or ragged ones are weak items.
#'
#' With correlated factors the bars of an item do not add up to its explained
#' variance exactly (the factors overlap), so the unique variance is taken
#' from the model, not from the bars.
#'
#' History: first drawn (2026-11-16) as the "grill" / "loom" plot with the
#' unique variance growing leftwards from the first line; the author found
#' that confusing and asked for a plain extra column.
#' @param L Data frame: item, own (the item's hypothesised factor), factor,
#'   and the standardized loading on `factor` in column `value`.
#' @param unique Named numeric vector: unique variance (standardized) per item.
#' @keywords internal
#' @noRd
psw_brush <- function(L, unique, value = "std", fit_text = NULL, fit_pos = c("none", "below", "above")) {
  fn <- unique(L$factor)
  m <- length(fn)
  items <- unique(L[c("item", "own")])
  items$y <- -(seq_len(nrow(items)) + (match(items$own, fn) - 1L) * 0.8)
  L$y <- items$y[match(L$item, items$item)]
  L$x0 <- (match(L$factor, fn) - 1L) * 1.08
  L$v <- L[[value]]
  L <- L[is.finite(L$v), , drop = FALSE]
  L$len <- pmin(L$v^2, 1)
  L$target <- L$factor == L$own
  items$u <- pmin(pmax(unname(unique[items$item]), 0), 1)
  ux <- m * 1.08 + 0.12                       # the unique-variance column, set a little apart
  ink <- "grey20"
  wl <- (max(nchar(psw_short(items$item, 22L))) * psw_cb$char + 14) / psw_cb$px_x
  handles <- data.frame(x = (seq_len(m) - 1L) * 1.08, lab = psw_short(fn, 16L))
  ytop <- max(items$y); ybot <- min(items$y)
  width <- round((ux + 1.15 + wl) * 150 + 30)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = items, ggplot2::aes(x = 0, xend = ux + 1, y = .data$y, yend = .data$y), colour = "grey92", linewidth = 0.3) +
    ggplot2::geom_segment(data = L, ggplot2::aes(x = .data$x0, xend = .data$x0 + .data$len, y = .data$y, yend = .data$y,
                                                 colour = .data$target), linewidth = 2.6) +
    ggplot2::geom_text(data = L[abs(L$v) >= 0.095, , drop = FALSE],
                       ggplot2::aes(x = .data$x0 + .data$len + 0.02, y = .data$y, label = psw_fmt_loading(.data$v)),
                       hjust = 0, size = 2.9, colour = ink) +
    ggplot2::geom_segment(data = items, ggplot2::aes(x = ux, xend = ux + .data$u, y = .data$y, yend = .data$y), colour = "grey80", linewidth = 2.6) +
    ggplot2::geom_text(data = items, ggplot2::aes(x = ux + .data$u + 0.02, y = .data$y, label = psw_fmt_loading(.data$u)),
                       hjust = 0, size = 2.9, colour = "grey45") +
    ggplot2::geom_text(data = items, ggplot2::aes(x = -0.05, y = .data$y, label = psw_short(.data$item, 22L)), hjust = 1, size = 3.4, colour = ink) +
    ggplot2::geom_segment(data = handles, ggplot2::aes(x = .data$x, xend = .data$x, y = ytop + 0.6, yend = ybot - 0.6),
                          colour = ink, linewidth = 2.2, lineend = "square") +
    ggplot2::annotate("segment", x = ux, xend = ux, y = ytop + 0.6, yend = ybot - 0.6, colour = "grey60", linewidth = 0.8) +
    ggplot2::geom_text(data = handles, ggplot2::aes(x = .data$x + 0.04, y = ytop + 1.1, label = .data$lab), hjust = 0,
                       fontface = "bold", size = 3.7, colour = ink) +
    ggplot2::annotate("text", x = ux + 0.04, y = ytop + 1.1, label = "Unique", hjust = 0, size = 3.6, colour = "grey45") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = ink, `FALSE` = "grey55"), guide = "none") +
    ggplot2::scale_x_continuous(breaks = c(handles$x, handles$x + 0.5, ux, ux + 0.5),
                                labels = c(rep("0", m), rep(".50", m), "0", ".50")) +
    ggplot2::coord_cartesian(xlim = c(-0.07 - wl, ux + 1.08), ylim = c(ybot - 0.8, ytop + 1.6), expand = FALSE, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    psw_theme_bars()
  attr(p, "psw_size") <- c(width = width, height = round((ytop - ybot + 2.4) * 22 + 65))
  psw_attach_fit(p, fit_text, fit_pos)
}

# ---- plot methods ------------------------------------------------------------

#' Figures for psyWiz results
#'
#' All figures have a fixed layout and are drawn in black, white and grey.
#'
#' * `"comb"`: one comb per factor, its name written along the spine. Each
#'   item is a tooth whose dark part is its squared loading (the variance it
#'   shares with the factor) and whose pale remainder is its unique variance,
#'   so every length is on the same 0 to 1 scale. Factor correlations are
#'   arcs joining the spines (thicker = stronger, up to four factors),
#'   residual correlations in the model are brackets at the tooth tips, and
#'   in a second-order model each factor has a top tooth showing how much of
#'   it the higher factor explains.
#' * `"path"`: factors on the left, items on the right, one arrow per
#'   loading; factor correlations as arcs (up to four factors).
#' * `"scalogram"` (`psw_scale`): respondents sorted by total score against
#'   the items; a Guttman scale shows as a clean staircase.
#' * `"antenna"`, `"facing"`, `"path"` (`psw_bifactor`): general factor on
#'   one side, specific factors on the other. In the antenna plot both grow
#'   from a mast in the middle; in the facing combs they grow towards each
#'   other and the gap between them is the item's unique variance.
#' * `"hasse"` (`psw_order`): the order found among items, the most common
#'   at the bottom; a line upwards means "implies".
#' * `"brush"` (`psw_esem`): every item on every factor, for models with
#'   cross-loadings. A clean factor looks like a toothbrush.
#'
#' There is one panel per group (`psw_invariance`) or per estimator
#' (`psw_multi`). The layout assumes a plotting area of roughly the size
#' given in the `psw_size` attribute of the returned plot (pixels).
#'
#' @param x A psyWiz result.
#' @param type Which figure, see above.
#' @param step For `psw_invariance`: which model of the ladder to draw, e.g.
#'   `"configural"` (default, every group has its own estimates), `"metric"`
#'   or `"scalar"`.
#' @param thick Comb plot: make the thickness of a tooth proportional to the
#'   observed variance of the item?
#' @param numbers,width Path diagram: print the loadings on the arrows; make
#'   the arrow width proportional to the loading.
#' @param order For `type = "scalogram"`: items in the order `"entered"` or
#'   sorted by `"difficulty"` (easiest first).
#' @param ... Ignored.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' plot(factorWiz(wizdata, f))
#' }
#' @name plot.psyWiz
NULL

#' @rdname plot.psyWiz
#' @export
plot.psw_cfa <- function(x, type = c("comb", "path"), thick = FALSE, numbers = TRUE, width = TRUE,
                         fit = c("none", "below", "above"), ...) {
  fit <- match.arg(fit)
  f <- psw_cfa_figure_data(x)
  fit_text <- if (fit != "none") psw_format_fit_line(x$fit, x$info$estimator) else NULL
  psw_figure(f$loadings, f$cors, "group", match.arg(type),
             opts = list(thick = thick, resid = f$resid, higher = x$higher, numbers = numbers, width = width,
                         fit_text = fit_text, fit_pos = fit))
}

# A single-group CFA drawn as one panel titled with the sample size.
psw_cfa_figure_data <- function(x) {
  panel <- sprintf("N = %d", x$info$n)
  l <- x$loadings
  l$group <- panel
  cc <- x$cors
  if (!is.null(cc) && nrow(cc)) cc$group <- panel else cc <- NULL
  rs <- x$residuals
  if (!is.null(rs) && nrow(rs)) rs$group <- panel else rs <- NULL
  list(loadings = l, cors = cc, resid = rs)
}

#' @rdname plot.psyWiz
#' @export
plot.psw_invariance <- function(x, type = c("path", "comb", "means"), step = "configural", thick = FALSE,
                                numbers = TRUE, width = TRUE, fit = c("none", "below", "above"), ...) {
  type <- match.arg(type)
  fit <- match.arg(fit)
  if (type == "means") {
    m <- if (!is.null(x$means)) x$means else invarianceMeans(x, ...)
    return(psw_plot_means(m, ...))
  }
  l <- x$loadings[x$loadings$step == step, , drop = FALSE]
  if (is.null(l) || nrow(l) == 0L)
    stop("No estimates for step '", step, "'. Available: ",
         paste(unique(x$loadings$step), collapse = ", "), call. = FALSE)
  cc <- if (is.null(x$cors)) NULL else x$cors[x$cors$step == step, , drop = FALSE]
  fit_row <- if (!is.null(x$fit)) x$fit[x$fit$step == step, , drop = FALSE] else NULL
  fit_text <- if (fit != "none") psw_format_fit_line(fit_row, x$info$estimator) else NULL
  psw_figure(l, cc, "group", type, caption = paste0(l$label[1L], " model."),
             opts = list(thick = thick, numbers = numbers, width = width, fit_text = fit_text, fit_pos = fit))
}

#' @rdname plot.psyWiz
#' @export
plot.psw_means <- function(x, ...) {
  psw_plot_means(x, ...)
}

#' Plot latent mean comparisons across groups (internal)
#'
#' @param means An object of class `psw_means` from `invarianceMeans()`.
#' @param fit_text Character string with formatted fit indices.
#' @param fit_pos Position of fit line (`"none"`, `"below"`, `"above"`).
#' @param ... Ignored.
#' @return A ggplot object with attribute `psw_size`.
#' @keywords internal
#' @noRd
psw_plot_means <- function(means, fit_text = NULL, fit_pos = c("none", "below", "above"), ...) {
  if (!isTRUE(attr(means, "supported"))) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "Latent means cannot be compared:\nscalar invariance is not supported across groups.",
                        size = 4, colour = "grey40", hjust = 0.5, vjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Latent Mean Comparisons")
    attr(p, "psw_size") <- c(width = 520, height = 240)
    return(p)
  }

  df <- as.data.frame(means)
  u_groups <- unique(df$group)
  ref <- attr(means, "reference")
  df$group <- factor(df$group, levels = rev(u_groups))

  alpha_val <- attr(means, "alpha")
  if (is.null(alpha_val)) alpha_val <- 0.05
  ci_pct <- round((1 - alpha_val) * 100)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mean, y = .data$group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
                           width = 0.22, linewidth = 0.65, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$is_ref, fill = .data$is_ref), size = 3.2, colour = "grey20") +
    ggplot2::scale_shape_manual(values = c("TRUE" = 22, "FALSE" = 21), guide = "none") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "grey80", "FALSE" = "grey20"), guide = "none") +
    ggplot2::facet_wrap(~ factor, scales = "free_x") +
    psw_theme_bars() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9.5, colour = "grey20", hjust = 1),
                   axis.ticks.y = ggplot2::element_line(colour = "grey80"),
                   plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey40", hjust = 0)) +
    ggplot2::labs(title = sprintf("Latent Mean Differences (%d%% CI)", ci_pct),
                  subtitle = sprintf("Reference group: %s (baseline = 0.00)%s", ref,
                                     if (identical(attr(means, "model_type"), "partial")) " [partial scalar]" else ""),
                  x = "Difference in latent mean", y = NULL)

  K <- length(unique(df$factor))
  G <- length(u_groups)
  w <- max(480, min(950, K * 260 + 60))
  h <- max(240, min(700, G * 36 + 120))
  attr(p, "psw_size") <- c(width = w, height = h)
  psw_attach_fit(p, fit_text, fit_pos)
}

#' @rdname plot.psyWiz
#' @export
plot.psw_multi <- function(x, type = c("path", "comb"), thick = FALSE, numbers = TRUE, width = TRUE, ...) {
  if (is.null(x$loadings)) stop("No model converged, nothing to draw.", call. = FALSE)
  psw_figure(x$loadings, x$cors, "estimator", match.arg(type), opts = list(thick = thick, numbers = numbers, width = width))
}

#' @rdname plot.psyWiz
#' @export
plot.psw_scale <- function(x, type = c("comb", "path", "scalogram"), thick = FALSE, numbers = TRUE, width = TRUE,
                           order = c("entered", "difficulty"), fit = c("none", "below", "above"), ...) {
  type <- match.arg(type)
  if (type == "scalogram") return(psw_scalogram_plot(x$guttman[[match.arg(order)]]$scalogram))
  if (is.null(x$loadings)) stop("The one-factor model did not converge, nothing to draw.", call. = FALSE)
  fit <- match.arg(fit)
  fit_text <- if (fit != "none" && !is.null(x$fit)) psw_format_fit_line(x$fit[1L, , drop = FALSE], x$info$estimator) else NULL
  l <- x$loadings
  l$group <- sprintf("N = %d", x$info$n)
  psw_figure(l, NULL, "group", type, opts = list(thick = thick, numbers = numbers, width = width,
                                                 fit_text = fit_text, fit_pos = fit))
}

#' @rdname plot.psyWiz
#' @export
plot.psw_bifactor <- function(x, type = c("antenna", "facing", "path"), numbers = TRUE, width = TRUE,
                              fit = c("none", "below", "above"), ...) {
  fit <- match.arg(fit)
  fit_text <- if (fit != "none") psw_format_fit_line(x$fit[1L, , drop = FALSE], x$info$estimator) else NULL
  switch(match.arg(type),
         antenna = psw_antenna(x$loadings, x$info$general, fit_text = fit_text, fit_pos = fit),
         facing = psw_facing(x$loadings, x$info$general, fit_text = fit_text, fit_pos = fit),
         path = psw_bifactor_path(x$loadings, x$info$general, numbers, width, fit_text = fit_text, fit_pos = fit))
}

#' @rdname plot.psyWiz
#' @export
plot.psw_esem <- function(x, type = "brush", fit = c("none", "below", "above"), ...) {
  if (all(is.na(x$loadings$esem))) stop("The ESEM did not converge, nothing to draw.", call. = FALSE)
  fit <- match.arg(fit)
  fit_row <- if (!is.null(x$fit)) x$fit[grepl("ESEM", x$fit$label), , drop = FALSE] else NULL
  if (is.null(fit_row) || nrow(fit_row) == 0L) fit_row <- x$fit[2L, , drop = FALSE]
  fit_text <- if (fit != "none") psw_format_fit_line(fit_row, x$info$estimator) else NULL
  psw_brush(x$loadings, x$uniqueness, value = "esem", fit_text = fit_text, fit_pos = fit)
}

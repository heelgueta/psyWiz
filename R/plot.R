# plot.R -- model diagrams.
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
#' @param highlight Colour weak loadings (below 0.40)?
#' @return A ggplot object.
#' @keywords internal
#' @noRd
psw_diagram <- function(loadings, cors = NULL, panel_col = "group", highlight = FALSE, caption = NULL) {
  lay <- psw_dg_layout(loadings, cors, panel_col)
  it <- lay$items
  panel <- function(d) factor(d[[panel_col]], levels = lay$panels)
  ink <- "grey25"
  weak_col <- "#B2182B"

  m <- match(paste(loadings$factor, loadings$item), paste(it$factor, it$item))
  arrows <- data.frame(panel = panel(loadings), x = 0.03, xend = 0.97,
                       y = unname(lay$fy[loadings$factor]), yend = it$y[m],
                       std = loadings$std, stringsAsFactors = FALSE)
  arrows <- arrows[is.finite(arrows$std), , drop = FALSE]
  arrows$lw <- 0.25 + 1.5 * pmin(abs(arrows$std), 1)
  arrows$col <- if (highlight) ifelse(abs(arrows$std) < psw_dg$weak, weak_col, ink) else ink
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
      lineend = "round", arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches"), type = "closed")) +
    ggplot2::geom_label(
      data = arrows, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$lab, colour = .data$col),
      size = 3.1, fill = "white", label.size = 0, label.padding = ggplot2::unit(0.08, "lines")) +
    ggplot2::geom_label(
      data = inodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      hjust = 0, size = 3.5, colour = ink, fill = "white", label.size = 0.3,
      label.r = ggplot2::unit(0, "lines"), label.padding = ggplot2::unit(0.3, "lines")) +
    ggplot2::geom_label(
      data = fnodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      hjust = 1, size = 3.7, fontface = "bold", colour = ink, fill = "grey93", label.size = 0.3,
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
  note <- "Standardized loadings; line width is proportional to the loading."
  if (highlight) note <- paste(note, sprintf("Loadings below %.2f in red.", psw_dg$weak))
  if (!lay$show_cors && !is.null(cors) && nrow(cors) > 0L)
    note <- paste(note, "Factor correlations are not drawn for more than four factors.")
  p +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_wrap(~panel, ncol = ncol) +
    ggplot2::coord_cartesian(xlim = lay$xlim, ylim = lay$ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(caption = if (is.null(caption)) note else paste(caption, note)) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 12,
                                                      margin = ggplot2::margin(4, 0, 6, 0)),
                   plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0,
                                                        margin = ggplot2::margin(8, 0, 0, 0)),
                   panel.spacing = ggplot2::unit(14, "pt"),
                   plot.margin = ggplot2::margin(8, 10, 8, 10),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

#' Model diagrams for psyWiz results
#'
#' Draws the fitted model with a fixed layout: factors on the left, items on
#' the right, one arrow per loading whose width is proportional to the
#' standardised loading, and factor correlations as arcs (for up to four
#' factors). There is one panel per group (`psw_invariance`) or per estimator
#' (`psw_multi`), so the panels can be compared directly.
#'
#' The layout assumes a plotting area of roughly the size given in the
#' `psw_size` attribute of the returned plot (pixels); much smaller devices
#' will crowd the labels.
#'
#' @param x A `psw_invariance` or `psw_multi` object.
#' @param step For `psw_invariance`: which model of the ladder to draw, e.g.
#'   `"configural"` (default, every group has its own estimates), `"metric"`
#'   or `"scalar"`.
#' @param highlight Colour loadings below 0.40 in red? Off by default.
#' @param ... Ignored.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' plot(cfaWiz(wizdata, f, estimators = c("MLR", "WLSMV")))
#' }
#' @name plot.psyWiz
NULL

#' @rdname plot.psyWiz
#' @export
plot.psw_invariance <- function(x, step = "configural", highlight = FALSE, ...) {
  l <- x$loadings[x$loadings$step == step, , drop = FALSE]
  if (is.null(l) || nrow(l) == 0L)
    stop("No estimates for step '", step, "'. Available: ",
         paste(unique(x$loadings$step), collapse = ", "), call. = FALSE)
  cc <- if (is.null(x$cors)) NULL else x$cors[x$cors$step == step, , drop = FALSE]
  p <- psw_diagram(l, cc, panel_col = "group", highlight = highlight,
                   caption = paste0(l$label[1L], " model."))
  attr(p, "psw_size") <- psw_diagram_size(l, cc, "group")
  p
}

#' @rdname plot.psyWiz
#' @export
plot.psw_multi <- function(x, highlight = FALSE, ...) {
  if (is.null(x$loadings)) stop("No model converged, nothing to draw.", call. = FALSE)
  p <- psw_diagram(x$loadings, x$cors, panel_col = "estimator", highlight = highlight)
  attr(p, "psw_size") <- psw_diagram_size(x$loadings, x$cors, "estimator")
  p
}

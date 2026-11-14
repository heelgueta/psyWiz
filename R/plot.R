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

psw_dg_layout <- function(loadings, cors, panel_col, tiles = FALSE) {
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
  # a variance tile is 0.8 rows high and (about) square on screen
  tw <- if (tiles) 0.8 * psw_dg$px_y / psw_dg$px_x + 0.04 else 0
  wi <- wi + tw
  list(items = items, fnames = fnames, fy = fy, flab = flab, ilab = ilab, wf = wf, wi = wi, tw = tw,
       show_cors = show_cors, arc = arc,
       xlim = c(-(wf + arc + if (show_cors) 0.52 else 0.02), 1 + wi + 0.04),
       ylim = c(min(items$y) - 0.7, -0.3),
       panels = unique(loadings[[panel_col]]))
}

# Image size in pixels that matches the layout assumptions above.
psw_diagram_size <- function(loadings, cors = NULL, panel_col = "group", tiles = FALSE) {
  lay <- psw_dg_layout(loadings, cors, panel_col, tiles)
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
psw_diagram <- function(loadings, cors = NULL, panel_col = "group", highlight = FALSE, caption = NULL,
                        tiles = FALSE) {
  lay <- psw_dg_layout(loadings, cors, panel_col, tiles)
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
      size = 3.1, fill = "white", label.size = 0, label.padding = ggplot2::unit(0.08, "lines"))
  if (tiles) {
    # variance tiles: the square is the item's variance, the filled square inside
    # it (same corner, area proportional to R-squared) is the part the factor explains
    side <- lay$tw - 0.04
    tl <- data.frame(panel = arrows$panel, y = arrows$yend, s = sqrt(pmin(arrows$std^2, 1)), col = arrows$col)
    p <- p +
      ggplot2::geom_rect(data = inodes, ggplot2::aes(xmin = 1, xmax = 1 + side, ymin = .data$y - 0.4,
                                                     ymax = .data$y + 0.4),
                         fill = "white", colour = ink, linewidth = 0.3) +
      ggplot2::geom_rect(data = tl, ggplot2::aes(xmin = 1, xmax = 1 + side * .data$s, ymin = .data$y - 0.4,
                                                 ymax = .data$y - 0.4 + 0.8 * .data$s, fill = .data$col)) +
      ggplot2::scale_fill_identity() +
      ggplot2::geom_text(data = inodes, ggplot2::aes(x = 1 + lay$tw, y = .data$y, label = .data$lab),
                         hjust = 0, size = 3.5, colour = ink)
  } else {
    p <- p + ggplot2::geom_label(
      data = inodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      hjust = 0, size = 3.5, colour = ink, fill = "white", label.size = 0.3,
      label.r = ggplot2::unit(0, "lines"), label.padding = ggplot2::unit(0.3, "lines"))
  }
  p <- p + ggplot2::geom_label(
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
  if (tiles) note <- paste(note, "Each square is an item's variance; the filled part is the share explained by the factor.")
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

# ---- comb plot ---------------------------------------------------------------
#
# One comb per factor. The spine is the factor; the top tooth has length 1,
# the factor's whole variance. Every item is a tooth whose dark part has
# length loading^2, the share of the item's variance it has in common with the
# factor, on a pale tooth of length 1 (so the pale remainder is the item's
# unique variance). Below the top tooth, one thin tooth per other factor shows
# the variance the two factors share (r^2). Everything in the figure is a
# length on the same 0-1 scale of shared variance, which is what makes combs
# from different groups, estimators or scales directly comparable.

psw_cb <- list(px_x = 210, px_y = 22, char = 7.2)

psw_comb_layout <- function(loadings, cors, panel_col) {
  items <- unique(loadings[c("factor", "item")])
  fnames <- unique(items$factor)
  show_cors <- !is.null(cors) && nrow(cors) > 0L && length(fnames) <= psw_dg$max_cor_factors
  nsub <- if (show_cors) length(fnames) - 1L else 0L
  y <- 0
  top <- numeric(length(fnames))
  sub <- list()
  items$y <- NA_real_
  for (f in seq_along(fnames)) {
    top[f] <- y
    for (g in setdiff(seq_along(fnames), f)[seq_len(nsub)]) {
      y <- y - 0.62
      sub[[length(sub) + 1L]] <- data.frame(factor = fnames[f], other = fnames[g], y = y, stringsAsFactors = FALSE)
    }
    y <- y - 0.35
    idx <- which(items$factor == fnames[f])
    items$y[idx] <- y - seq_along(idx)
    y <- min(items$y[idx]) - 1.6
  }
  ilab <- psw_short(items$item, 22L)
  flab <- psw_short(fnames, 22L)
  left <- (max(nchar(c(ilab, flab))) * psw_cb$char + 14) / psw_cb$px_x
  right <- if (show_cors) 0.38 else 0.17
  list(items = items, ilab = ilab, fnames = fnames, flab = flab, top = stats::setNames(top, fnames),
       sub = if (length(sub)) do.call(rbind, sub) else NULL, show_cors = show_cors,
       xlim = c(-left, 1 + right), ylim = c(min(items$y) - 0.8, 0.7),
       panels = unique(loadings[[panel_col]]))
}

psw_comb_size <- function(loadings, cors = NULL, panel_col = "group") {
  lay <- psw_comb_layout(loadings, cors, panel_col)
  np <- length(lay$panels)
  ncol <- if (np <= 3L) np else 3L
  nrow <- ceiling(np / ncol)
  c(width = max(520, round(ncol * diff(lay$xlim) * psw_cb$px_x + 30)),
    height = round(nrow * (diff(lay$ylim) * psw_cb$px_y + 28) + 80))
}

#' Comb plot of a simple-structure factor model (internal)
#' @inheritParams psw_diagram
#' @keywords internal
#' @noRd
psw_comb <- function(loadings, cors = NULL, panel_col = "group", highlight = FALSE, caption = NULL) {
  lay <- psw_comb_layout(loadings, cors, panel_col)
  it <- lay$items
  lv <- lay$panels
  ink <- "grey20"
  weak_col <- "#B2182B"
  panel <- function(d) factor(d[[panel_col]], levels = lv)

  m <- match(paste(loadings$factor, loadings$item), paste(it$factor, it$item))
  teeth <- data.frame(panel = panel(loadings), y = it$y[m], std = loadings$std, stringsAsFactors = FALSE)
  teeth <- teeth[is.finite(teeth$std), , drop = FALSE]
  teeth$len <- pmin(teeth$std^2, 1)
  teeth$col <- if (highlight) ifelse(abs(teeth$std) < psw_dg$weak, weak_col, ink) else ink
  teeth$lab <- psw_fmt_loading(teeth$std)

  pg <- function(n) factor(rep(lv, each = n), levels = lv)
  nf <- length(lay$fnames)
  ghost <- data.frame(panel = pg(nrow(it)), y = rep(it$y, length(lv)))
  inames <- data.frame(panel = pg(nrow(it)), y = rep(it$y, length(lv)), lab = rep(lay$ilab, length(lv)))
  tops <- data.frame(panel = pg(nf), y = rep(unname(lay$top), length(lv)), lab = rep(lay$flab, length(lv)))
  bottom <- vapply(lay$fnames, function(f) min(it$y[it$factor == f]), numeric(1))
  spine <- data.frame(panel = pg(nf), y = rep(unname(lay$top), length(lv)), yend = rep(unname(bottom), length(lv)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ghost, ggplot2::aes(x = 0, xend = 1, y = .data$y, yend = .data$y),
                          colour = "grey88", linewidth = 2.6) +
    ggplot2::geom_segment(data = teeth, ggplot2::aes(x = 0, xend = .data$len, y = .data$y, yend = .data$y,
                                                     colour = .data$col), linewidth = 2.6) +
    ggplot2::geom_text(data = teeth, ggplot2::aes(x = .data$len + 0.02, y = .data$y, label = .data$lab,
                                                  colour = .data$col), hjust = 0, size = 3.2) +
    ggplot2::geom_text(data = inames, ggplot2::aes(x = -0.035, y = .data$y, label = .data$lab),
                       hjust = 1, size = 3.4, colour = ink) +
    ggplot2::geom_segment(data = tops, ggplot2::aes(x = 0, xend = 1, y = .data$y, yend = .data$y),
                          colour = ink, linewidth = 3.4) +
    ggplot2::geom_text(data = tops, ggplot2::aes(x = -0.035, y = .data$y, label = .data$lab),
                       hjust = 1, size = 3.7, fontface = "bold", colour = ink)

  if (lay$show_cors) {
    both <- rbind(data.frame(panel = panel(cors), factor = cors$factor1, other = cors$factor2, r = cors$r),
                  data.frame(panel = panel(cors), factor = cors$factor2, other = cors$factor1, r = cors$r))
    both$y <- lay$sub$y[match(paste(both$factor, both$other), paste(lay$sub$factor, lay$sub$other))]
    both <- both[is.finite(both$y) & is.finite(both$r), , drop = FALSE]
    both$len <- pmin(both$r^2, 1)
    both$lab <- sprintf("%s  r = %s", psw_short(both$other, 14L), psw_fmt_loading(both$r))
    p <- p +
      ggplot2::geom_segment(data = both, ggplot2::aes(x = 0, xend = .data$len, y = .data$y, yend = .data$y),
                            colour = "grey50", linewidth = 1.5) +
      ggplot2::geom_text(data = both, ggplot2::aes(x = .data$len + 0.02, y = .data$y, label = .data$lab),
                         hjust = 0, size = 2.9, colour = "grey35")
  }

  ncol <- if (length(lv) <= 3L) length(lv) else 3L
  note <- paste0("Tooth length is shared variance. Top tooth: the factor (1). Items: squared standardized ",
                 "loading (the number is the loading); the pale part is unique variance.")
  if (lay$show_cors) note <- paste(note, "Thin teeth: variance shared with another factor (r squared).")
  if (highlight) note <- paste(note, sprintf("Loadings below %.2f in red.", psw_dg$weak))
  p +
    ggplot2::geom_segment(data = spine, ggplot2::aes(x = 0, xend = 0, y = .data$y, yend = .data$yend),
                          colour = ink, linewidth = 3.4, lineend = "square") +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::facet_wrap(~panel, ncol = ncol) +
    ggplot2::coord_cartesian(xlim = lay$xlim, ylim = lay$ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste(strwrap(if (is.null(caption)) note else paste(caption, note),
                                          width = max(85, round(ncol * diff(lay$xlim) * psw_cb$px_x / 5.6))),
                                  collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
                   axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 8.5, colour = "grey40"),
                   strip.text = ggplot2::element_text(face = "bold", size = 12),
                   plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0,
                                                        margin = ggplot2::margin(8, 0, 0, 0)),
                   panel.spacing = ggplot2::unit(16, "pt"), plot.margin = ggplot2::margin(8, 12, 8, 12),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# One entry point for the three figure types.
psw_figure <- function(loadings, cors, panel_col, type = "diagram", highlight = FALSE, caption = NULL) {
  if (type == "comb") {
    p <- psw_comb(loadings, cors, panel_col, highlight, caption)
    attr(p, "psw_size") <- psw_comb_size(loadings, cors, panel_col)
  } else {
    tiles <- type == "tiles"
    p <- psw_diagram(loadings, cors, panel_col, highlight, caption, tiles = tiles)
    attr(p, "psw_size") <- psw_diagram_size(loadings, cors, panel_col, tiles)
  }
  p
}

#' Figures for psyWiz results
#'
#' Three ways of drawing a fitted factor model, all with a fixed layout and
#' one panel per group (`psw_invariance`), per estimator (`psw_multi`) or for
#' the single scale (`psw_scale`), so panels can be compared directly.
#'
#' * `"diagram"`: factors on the left, items on the right, one arrow per
#'   loading whose width is proportional to the standardised loading; factor
#'   correlations as arcs (up to four factors).
#' * `"tiles"`: the same, with a square next to every item. The square is the
#'   item's variance and the filled square inside it (area proportional to
#'   R-squared) is the part explained by the factor.
#' * `"comb"`: one comb per factor. The top tooth is the factor's variance
#'   (length 1); each item is a tooth whose dark part is its squared loading,
#'   the variance it shares with the factor, and whose pale remainder is its
#'   unique variance. Thin teeth under the top tooth show the variance shared
#'   with each other factor (r squared). Every length is on the same 0 to 1
#'   scale of shared variance.
#'
#' The layout assumes a plotting area of roughly the size given in the
#' `psw_size` attribute of the returned plot (pixels); much smaller devices
#' will crowd the labels.
#'
#' @param x A `psw_invariance`, `psw_multi` or `psw_scale` object.
#' @param type `"diagram"`, `"tiles"` or `"comb"`.
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
#' plot(cfaWiz(wizdata, f, estimators = c("MLR", "WLSMV")), type = "comb")
#' }
#' @name plot.psyWiz
NULL

#' @rdname plot.psyWiz
#' @export
plot.psw_invariance <- function(x, type = c("diagram", "tiles", "comb"), step = "configural",
                                highlight = FALSE, ...) {
  type <- match.arg(type)
  l <- x$loadings[x$loadings$step == step, , drop = FALSE]
  if (is.null(l) || nrow(l) == 0L)
    stop("No estimates for step '", step, "'. Available: ",
         paste(unique(x$loadings$step), collapse = ", "), call. = FALSE)
  cc <- if (is.null(x$cors)) NULL else x$cors[x$cors$step == step, , drop = FALSE]
  psw_figure(l, cc, "group", type, highlight, caption = paste0(l$label[1L], " model."))
}

#' @rdname plot.psyWiz
#' @export
plot.psw_multi <- function(x, type = c("diagram", "tiles", "comb"), highlight = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(x$loadings)) stop("No model converged, nothing to draw.", call. = FALSE)
  psw_figure(x$loadings, x$cors, "estimator", type, highlight)
}

#' @rdname plot.psyWiz
#' @export
plot.psw_scale <- function(x, type = c("comb", "diagram", "tiles"), highlight = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(x$loadings)) stop("The congeneric model did not converge, nothing to draw.", call. = FALSE)
  l <- x$loadings
  l$group <- sprintf("%s (N = %d)", x$info$name, x$info$n)
  psw_figure(l, NULL, "group", type, highlight, caption = "Congeneric model.")
}

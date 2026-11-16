# engine-order.R -- exploratory order among binary items: which items form a
# cumulative line, which are interchangeable, which branch, which stand apart.
#
# Factor analysis asks what items have in common. Some item sets are better
# described by ORDER: whoever endorses a hard item also endorses the easier
# ones (Guttman, 1950). Real data are rarely one clean line; there can be
# several items on the same rung, branches (two items that both build on a
# third but not on each other), and items that belong to no order at all.
# This file looks for that structure without being told what it is.
#
# Three views of the same data, all computed here without other packages:
#
#   1. Scales (Mokken, 1971; Sijtsma & Molenaar, 2002). Items are grouped
#      bottom-up into cumulative scales by Loevinger's H, and items that fit
#      nowhere are listed as unscalable. A simplified version of Mokken's
#      automated item selection: no significance tests.
#   2. Implications (the idea behind item tree analysis: van Leeuwe, 1974;
#      Schrepp, 1999, 2003; Sargin & Unlu, 2009). "i implies j" is accepted
#      when few of those who endorse i fail to endorse j (at most `tolerance`)
#      AND that is far fewer than chance would give (pairwise H of at least
#      `h.min`). The second condition matters: without it, any item implies
#      any very common item. Items that imply each other share a rung. The
#      accepted implications form a partial order, drawn as a Hasse diagram.
#   3. Response patterns. How many respondents fit a single line, and how
#      many fit the partial order found in 2? The gap is the evidence that
#      the structure is more than a line.
#
# Not done, on purpose: POSAC (Shye, 1985) places response patterns in two
# dimensions by an iterative fit; the Hasse diagram answers the same question
# ("a line, a tree, or unrelated?") at the level of items and is far easier
# to read. The packages mokken, DAKS and POSAC cannot be used from here: a
# jamovi module has to bundle its dependencies, and none of them ships with
# jamovi.

# 0/1 matrix from items of any type. Two-valued items: the higher value is 1.
# Items with more values: 1 when the response is above the midpoint of the
# observed range. Constant items are dropped.
psw_dichotomise <- function(data, items) {
  missing_items <- setdiff(items, names(data))
  if (length(missing_items)) stop("Not in the data: ", paste(missing_items, collapse = ", "), call. = FALSE)
  X <- vapply(items, function(v) psw_as_numeric(data[[v]]), numeric(nrow(data)))
  X <- X[stats::complete.cases(X), , drop = FALSE]
  notes <- character()
  nval <- apply(X, 2L, function(v) length(unique(v)))
  if (any(nval > 2L))
    notes <- c(notes, sprintf("Items with more than two values were split at the midpoint of their range (1 = above it): %s.",
                              paste(items[nval > 2L], collapse = ", ")))
  B <- apply(X, 2L, function(v) as.numeric(v > (min(v) + max(v)) / 2))
  dim(B) <- dim(X)
  colnames(B) <- items
  const <- apply(B, 2L, function(v) length(unique(v)) < 2L)
  if (any(const))
    notes <- c(notes, sprintf("Left out because everybody gave the same answer: %s.", paste(items[const], collapse = ", ")))
  list(X = B[, !const, drop = FALSE], notes = notes, n.dropped = nrow(data) - nrow(X))
}

# Bottom-up grouping of items into cumulative scales (after Mokken's automated
# item selection procedure). A scale starts with the pair with the highest
# pairwise H; then the item with the highest H with respect to the scale is
# added, as long as that H is at least `c` and the item correlates positively
# with every item already in. Repeats with the items left over.
# Returns an integer vector: scale number per item, 0 = unscalable.
psw_aisp <- function(X, c = 0.3) {
  k <- ncol(X)
  S <- stats::cov(X)
  Smax <- stats::cov(apply(X, 2L, sort))
  Hij <- S / Smax
  diag(Hij) <- NA_real_
  member <- integer(k)
  left <- seq_len(k)
  s <- 0L
  while (length(left) >= 2L) {
    sub <- Hij[left, left, drop = FALSE]
    if (!any(is.finite(sub)) || max(sub, na.rm = TRUE) < c) break
    start <- left[arrayInd(which.max(sub), dim(sub))[1L, ]]
    scale <- start
    repeat {
      cand <- setdiff(left, scale)
      if (!length(cand)) break
      h <- vapply(cand, function(j)
        if (all(Hij[j, scale] > 0)) sum(S[j, scale]) / sum(Smax[j, scale]) else -Inf, numeric(1))
      if (max(h) < c) break
      scale <- c(scale, cand[which.max(h)])
    }
    s <- s + 1L
    member[scale] <- s
    left <- setdiff(left, scale)
  }
  member
}

#' Exploratory order among items
#'
#' For binary items (others are split at the midpoint of their range): finds
#' which items form a cumulative line, which share a rung, which branch and
#' which stand apart, without a hypothesis. See the comments at the top of
#' `R/engine-order.R` and `docs/METHODS.md` for the three views it combines.
#'
#' @param data A data frame.
#' @param items Character vector of item names (at least three).
#' @param tolerance Largest share of counterexamples for an implication:
#'   "i implies j" is accepted when at most this share of those who endorse
#'   i do not endorse j. Default 0.10.
#' @param h.min Smallest pairwise H (Loevinger) for an implication, and
#'   smallest item H for joining a scale. Default 0.30.
#' @return An object of class `psw_order`: a list with data frames `items`
#'   (endorsement, level in the order, items on the same rung, scale, H,
#'   role), `implications` (the direct ones), `all.implications`, `scales`,
#'   `patterns`, `fit` (how many respondents fit a line and the order), the
#'   `hasse` diagram data, a `reading` in words, plus `notes` and `info`.
#' @references
#' Guttman, L. (1950). The basis for scalogram analysis. In S. A. Stouffer et
#' al., *Measurement and prediction* (pp. 60-90). Princeton University Press.
#'
#' Mokken, R. J. (1971). *A theory and procedure of scale analysis*. Mouton.
#'
#' Sijtsma, K., & Molenaar, I. W. (2002). *Introduction to nonparametric item
#' response theory*. Sage.
#'
#' Schrepp, M. (2003). A method for the analysis of hierarchical dependencies
#' between items of a questionnaire. *Methods of Psychological Research
#' Online, 8*(1), 43-79.
#'
#' Sargin, A., & Unlu, A. (2009). Inductive item tree analysis: Corrections,
#' improvements, and comparisons. *Mathematical Social Sciences, 58*(3),
#' 376-392.
#'
#' Shye, S. (1985). *Multiple scaling: The theory and application of partial
#' order scalogram analysis*. North-Holland.
#' @examples
#' data(wizlab)
#' orderWiz(wizlab, grep("^debt_", names(wizlab), value = TRUE))
#' @export
orderWiz <- function(data, items, tolerance = 0.10, h.min = 0.30) {
  items <- unique(as.character(items))
  if (length(items) < 3L) stop("Choose at least three items.", call. = FALSE)
  dz <- psw_dichotomise(as.data.frame(data), items)
  X <- dz$X
  notes <- dz$notes
  if (dz$n.dropped > 0L) notes <- c(notes, sprintf("%d cases with missing values were left out.", dz$n.dropped))
  k <- ncol(X)
  n <- nrow(X)
  if (k < 3L) stop("Fewer than three items vary; there is no order to look for.", call. = FALSE)
  if (n < 20L) stop("Fewer than 20 complete cases.", call. = FALSE)
  nm <- colnames(X)
  p <- colMeans(X)
  mk <- psw_mokken_h(X)

  # ---- 1. scales ----
  member <- psw_aisp(X, h.min)
  scales <- do.call(rbind, lapply(sort(unique(member[member > 0L])), function(s) {
    idx <- which(member == s)
    data.frame(scale = s, k = length(idx), H = psw_mokken_h(X[, idx, drop = FALSE])$H,
               items = paste(nm[idx][order(-p[idx])], collapse = ", "), stringsAsFactors = FALSE)
  }))

  # ---- 2. implications: V[i, j] = share of those endorsing i who do not endorse j ----
  V <- vapply(seq_len(k), function(j) vapply(seq_len(k), function(i) sum(X[, i] == 1 & X[, j] == 0) / sum(X[, i]), numeric(1)), numeric(k))
  holds <- V <= tolerance & mk$Hij >= h.min
  holds[is.na(holds)] <- FALSE
  diag(holds) <- FALSE
  # rungs: items that imply each other
  mutual <- holds & t(holds)
  rung <- seq_len(k)
  for (i in seq_len(k)) for (j in seq_len(k)) if (mutual[i, j]) rung[rung == rung[j]] <- rung[i]
  rung <- match(rung, unique(rung))
  nr <- max(rung)
  rp <- vapply(seq_len(nr), function(r) mean(p[rung == r]), numeric(1))
  # rung a implies rung b when every item of a implies every item of b; only towards more common rungs
  R <- matrix(FALSE, nr, nr)
  for (a in seq_len(nr)) for (b in seq_len(nr))
    if (a != b && rp[a] < rp[b]) R[a, b] <- all(holds[rung == a, rung == b, drop = FALSE])
  # level = length of the longest chain of prerequisites below a rung
  level <- rep(1L, nr)
  for (a in order(-rp)) if (any(R[a, ])) level[a] <- 1L + max(level[R[a, ]])
  connected <- rowSums(R) + colSums(R) > 0 | tabulate(rung, nr) > 1L
  # direct implications only (drop a -> c when a -> b -> c)
  cover <- R
  for (a in seq_len(nr)) for (b in which(R[a, ])) if (any(R[a, ] & R[, b])) cover[a, b] <- FALSE

  rung_label <- vapply(seq_len(nr), function(r) paste(nm[rung == r][order(-p[rung == r])], collapse = " = "), "")
  pair_stats <- function(M) {
    idx <- which(M, arr.ind = TRUE)
    if (!nrow(idx)) return(NULL)
    out <- do.call(rbind, lapply(seq_len(nrow(idx)), function(r) {
      a <- which(rung == idx[r, 1L]); b <- which(rung == idx[r, 2L])
      data.frame(from = rung_label[idx[r, 1L]], to = rung_label[idx[r, 2L]],
                 exceptions = max(V[a, b]), H = min(mk$Hij[a, b]), stringsAsFactors = FALSE)
    }))
    out[order(-rp[idx[, 2L]], -rp[idx[, 1L]]), , drop = FALSE]
  }

  role <- ifelse(!connected[rung], if (any(connected)) "apart: no order with the other items" else "apart",
          ifelse(tabulate(rung, nr)[rung] > 1L, "shares a rung", "in the order"))
  items_df <- data.frame(
    item = nm, endorsed = unname(p), level = ifelse(connected[rung], level[rung], NA_integer_),
    same.as = vapply(seq_len(k), function(i) paste(setdiff(nm[rung == rung[i]], nm[i]), collapse = ", "), ""),
    scale = member, H = unname(mk$Hi), role = role, stringsAsFactors = FALSE)
  items_df <- items_df[order(is.na(items_df$level), items_df$level, -items_df$endorsed), , drop = FALSE]
  rownames(items_df) <- NULL

  # ---- 3. response patterns ----
  ord <- order(-p)
  Xo <- X[, ord, drop = FALSE]
  line_ok <- apply(Xo, 1L, function(v) !is.unsorted(rev(v)))          # 1..1 0..0 from the most common item on
  viol <- which(holds, arr.ind = TRUE)                                  # i endorsed without j
  order_ok <- if (nrow(viol)) rowSums(X[, viol[, 1L], drop = FALSE] == 1 & X[, viol[, 2L], drop = FALSE] == 0) == 0 else rep(TRUE, n)
  key <- apply(Xo, 1L, paste, collapse = "")
  tab <- sort(table(key), decreasing = TRUE)
  first <- match(names(tab), key)
  patterns <- data.frame(pattern = names(tab), n = as.integer(tab), share = as.numeric(tab) / n,
                         line = line_ok[first], order = order_ok[first], stringsAsFactors = FALSE)
  gt <- psw_guttman(X, "difficulty")
  fit <- data.frame(
    model = c("One line through all items", "The order found"),
    fitting = c(mean(line_ok), if (any(R) || any(mutual)) mean(order_ok) else NA_real_),
    detail = c(sprintf("coefficient of reproducibility %.3f, scalability %.3f", gt$cr, gt$cs),
               sprintf("%d implication(s) among %d rung(s); %d item(s) apart", sum(R), sum(connected), sum(!connected[rung]))),
    stringsAsFactors = FALSE)

  # ---- reading ----
  branching <- any(vapply(which(connected), function(a) {
    b <- setdiff(which(connected), a)
    any(!R[a, b] & !R[b, a])                       # a connected rung it is not comparable with
  }, logical(1)))
  apart <- nm[!connected[rung]]
  reading <- c(
    if (!any(connected)) "No implication between any two items was found: these items do not form an order."
    else if (!branching && !length(apart)) "All items fall on one line: a cumulative (Guttman-type) scale."
    else if (!branching) "The items in the order fall on one line."
    else "The order branches: some items build on the same items but not on each other, so no single line describes them.",
    if (any(tabulate(rung, nr) > 1L))
      sprintf("On the same rung (interchangeable): %s.", paste(rung_label[tabulate(rung, nr) > 1L], collapse = "; ")),
    if (length(apart) && any(connected)) sprintf("Apart from the order: %s.", paste(apart, collapse = ", ")),
    if (any(connected)) sprintf("%.0f%% of respondents fit one line through all items; %.0f%% fit the order found.",
                                100 * mean(line_ok), 100 * mean(order_ok))
    else sprintf("%.0f%% of respondents fit one line through all items.", 100 * mean(line_ok)))

  # ---- Hasse diagram: rungs as nodes, direct implications as edges ----
  nodes <- data.frame(id = seq_len(nr), label = rung_label, endorsed = rp, level = level, connected = connected,
                      stringsAsFactors = FALSE)
  e <- which(cover, arr.ind = TRUE)
  out <- list(
    items = items_df, implications = pair_stats(cover), all.implications = pair_stats(R),
    scales = scales, patterns = patterns, fit = fit, reading = reading,
    hasse = list(nodes = nodes, edges = data.frame(from = e[, 1L], to = e[, 2L])),
    notes = unique(notes),
    info = list(n = n, k = k, tolerance = tolerance, h.min = h.min, H = mk$H,
                order = nm[ord], unscalable = nm[member == 0L]))
  class(out) <- "psw_order"
  out
}

#' Hasse diagram of the order found by orderWiz() (internal)
#'
#' Rungs are boxes, the most common at the bottom. A line going up from A to
#' B means "whoever endorses B also endorses A" (direct implications only).
#' Boxes on the same height that are not joined are alternatives: a branch.
#' Items that are in no order with the others stand to the right.
#' @keywords internal
#' @noRd
psw_hasse <- function(h) {
  nd <- h$nodes
  ink <- "grey20"
  nd$lab <- sprintf("%s\n%.0f%%", gsub(" = ", " =\n", nd$label, fixed = TRUE), 100 * nd$endorsed)
  inn <- nd[nd$connected, , drop = FALSE]
  out <- nd[!nd$connected, , drop = FALSE]
  width <- function(lab) max(nchar(unlist(strsplit(lab, "\n")))) * 0.075 + 0.25
  nd$x <- NA_real_; nd$y <- NA_real_
  span <- 0
  for (lv in sort(unique(inn$level))) {
    idx <- which(nd$connected & nd$level == lv)
    idx <- idx[order(-nd$endorsed[idx])]
    w <- vapply(nd$lab[idx], width, numeric(1))
    centres <- cumsum(w) - w / 2 + 0.15 * (seq_along(idx) - 1)
    nd$x[idx] <- centres - (sum(w) + 0.15 * (length(idx) - 1)) / 2
    nd$y[idx] <- lv
    span <- max(span, sum(w) + 0.15 * (length(idx) - 1))
  }
  top <- max(c(1, inn$level))
  if (nrow(out)) {
    idx <- which(!nd$connected)
    idx <- idx[order(-nd$endorsed[idx])]
    wo <- max(vapply(nd$lab[idx], width, numeric(1)))
    nd$x[idx] <- span / 2 + 0.5 + wo / 2
    nd$y[idx] <- if (length(idx) == 1L) (1 + top) / 2 else seq(top, 1, length.out = length(idx))
    if (!nrow(inn)) nd$x[idx] <- 0
  } else wo <- 0
  ed <- h$edges
  seg <- data.frame(x = nd$x[ed$from], y = nd$y[ed$from], xend = nd$x[ed$to], yend = nd$y[ed$to])
  xl <- c(-span / 2 - 0.2, max(span / 2, if (nrow(out)) span / 2 + 0.5 + wo else 0) + 0.2)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                          colour = "grey45", linewidth = 0.7) +
    ggplot2::geom_label(data = nd[nd$connected, , drop = FALSE], ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                        size = 3.5, colour = ink, fill = "white", lineheight = 0.95,
                        label.padding = ggplot2::unit(0.35, "lines")) +
    ggplot2::geom_label(data = nd[!nd$connected, , drop = FALSE], ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                        size = 3.5, colour = "grey35", fill = "grey94", lineheight = 0.95,
                        label.padding = ggplot2::unit(0.35, "lines"))
  if (nrow(out) && nrow(inn))
    p <- p + ggplot2::annotate("text", x = nd$x[!nd$connected][1L], y = top + 0.55, label = "apart", size = 3.4,
                               colour = "grey45", fontface = "italic")
  w_px <- round(max(420, diff(xl) * 105 + 40))
  note <- paste0("Read upwards: whoever endorses an item also endorses the items it is joined to below. Items joined by = ",
                 "are interchangeable. Percentages: share endorsing. Items not joined, on any height, do not depend on each other.")
  p <- p +
    ggplot2::coord_cartesian(xlim = xl, ylim = c(0.5, top + 0.75), expand = FALSE, clip = "off") +
    ggplot2::labs(caption = paste(strwrap(note, floor((w_px - 30) / 6.3)), collapse = "\n")) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, colour = "grey40", hjust = 0, margin = ggplot2::margin(10, 0, 0, 0)),
                   plot.margin = ggplot2::margin(10, 12, 8, 12), plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  attr(p, "psw_size") <- c(width = w_px, height = round(top * 95 + 120))
  p
}

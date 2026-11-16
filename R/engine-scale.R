# engine-scale.R -- everything about ONE scale (one factor): reliability,
# the congeneric / tau-equivalent / parallel ladder, item statistics and
# Guttman / Mokken scalability.

# ---- classical statistics on a complete numeric item matrix ----------------

psw_alpha <- function(C) {
  k <- ncol(C)
  if (k < 2L) return(NA_real_)
  k / (k - 1) * (1 - sum(diag(C)) / sum(C))
}

# Guttman's lambda-2: never below alpha, and still a lower bound to reliability.
psw_lambda2 <- function(C) {
  k <- ncol(C)
  off <- C
  diag(off) <- 0
  (sum(off) + sqrt(k / (k - 1) * sum(off^2))) / sum(C)
}

#' Loevinger's H scalability coefficients (Mokken scale analysis), internal
#'
#' H compares the observed covariance of two items with the largest
#' covariance their marginal distributions allow (obtained when both items
#' order the respondents identically, i.e. a perfect Guttman scale):
#' H_ij = cov(X_i, X_j) / cov_max(X_i, X_j). For binary items this is 1 minus
#' the ratio of observed to expected Guttman errors; for polytomous items it
#' equals Molenaar's (1991) weighted-error definition. H = 1 is a perfect
#' (polytomous) Guttman scale, H = 0 is what unrelated items give.
#'
#' @param X Numeric matrix of complete item responses (higher = more).
#' @return list(H, Hi, Hij)
#' @keywords internal
#' @noRd
psw_mokken_h <- function(X) {
  X <- as.matrix(X)
  k <- ncol(X)
  S <- stats::cov(X)
  sorted <- apply(X, 2L, sort)
  Smax <- stats::cov(sorted)
  off <- row(S) != col(S)
  Hij <- S / Smax
  diag(Hij) <- NA_real_
  Hi <- vapply(seq_len(k), function(i) sum(S[i, -i]) / sum(Smax[i, -i]), numeric(1))
  list(H = sum(S[off]) / sum(Smax[off]), Hi = stats::setNames(Hi, colnames(X)), Hij = Hij)
}

#' Guttman scalogram statistics on item steps (internal)
#'
#' Each item with m categories is split into m - 1 binary "item steps"
#' (response >= category k). Steps are ordered from the most to the least
#' often passed; a respondent who passed s steps is expected to have passed
#' exactly the s easiest ones. Errors are counted with the Goodenough-Edwards
#' method (every cell that differs from that ideal pattern).
#'
#' CR  = coefficient of reproducibility = 1 - errors / (N * steps)
#' MMR = minimal marginal reproducibility: the CR obtainable from the
#'       marginals alone (mean over steps of max(p, 1 - p))
#' CS  = coefficient of scalability = (CR - MMR) / (1 - MMR)
#'
#' For binary items this is the classical Guttman analysis. For polytomous
#' items the steps of one item are nested by construction, which pushes CR
#' up; CS and Loevinger's H are the more informative numbers there.
#' @keywords internal
#' @noRd
psw_guttman <- function(X, order = c("difficulty", "entered")) {
  order <- match.arg(order)
  X <- as.matrix(X)
  steps <- do.call(cbind, lapply(seq_len(ncol(X)), function(j) {
    v <- sort(unique(X[, j]))[-1L]
    vapply(v, function(cut) as.numeric(X[, j] >= cut), numeric(nrow(X)))
  }))
  p <- colMeans(steps)
  # The order the scale is tested against. "difficulty": easiest step first, as
  # observed (the usual, exploratory scalogram). "entered": the order in which
  # the items were given, easiest first, i.e. a hypothesis. That is only
  # defined for binary items; steps of polytomous items are always sorted.
  binary <- ncol(steps) == ncol(X)
  ord <- if (order == "entered" && binary) seq_len(ncol(steps)) else order(-p)
  steps <- steps[, ord, drop = FALSE]
  p <- p[ord]
  s <- rowSums(steps)
  ideal <- outer(s, seq_len(ncol(steps)), ">=") * 1
  errors <- sum(steps != ideal)
  cr <- 1 - errors / length(steps)
  mmr <- mean(pmax(p, 1 - p))
  list(cr = cr, mmr = mmr, cs = (cr - mmr) / (1 - mmr), errors = errors,
       nsteps = ncol(steps), perfect = mean(rowSums(steps != ideal) == 0),
       order = if (binary) order else "difficulty")
}

#' Data for a scalogram figure (internal)
#'
#' Items in columns (in the entered order, or easiest first), respondents in
#' rows sorted by total score. Responses are rescaled to 0-1 per item. To keep
#' the figure light the sorted respondents are averaged in at most `bins`
#' consecutive groups, so a cell is the share of that group endorsing the item
#' (binary items) or its mean response (polytomous items).
#' @keywords internal
#' @noRd
psw_scalogram <- function(X, order = c("difficulty", "entered"), bins = 150L) {
  order <- match.arg(order)
  X <- as.matrix(X)
  lo <- apply(X, 2L, min)
  rng <- pmax(apply(X, 2L, max) - lo, 1e-12)
  Z <- sweep(sweep(X, 2L, lo), 2L, rng, `/`)
  io <- if (order == "entered") seq_len(ncol(Z)) else order(-colMeans(Z))
  Z <- Z[order(rowSums(Z)), io, drop = FALSE]
  nb <- min(nrow(Z), bins)
  grp <- ceiling(seq_len(nrow(Z)) * nb / nrow(Z))
  z <- rowsum(Z, grp) / as.numeric(table(grp))
  list(z = z, items = colnames(X)[io], binary = all(apply(X, 2L, function(v) length(unique(v)) <= 2L)),
       n = nrow(X), order = order)
}

# The scalability table and scalogram for one ordering of the items.
psw_scalability <- function(X, mk, order) {
  gt <- psw_guttman(X, order)
  hij <- mk$Hij[lower.tri(mk$Hij)]
  tab <- data.frame(
    statistic = c("Loevinger's H (scale)", "Smallest item H", "Smallest item-pair H",
                  "Coefficient of reproducibility (CR)", "Minimal marginal reproducibility (MMR)",
                  "Coefficient of scalability (CS)", "Respondents with a perfect Guttman pattern"),
    value = c(mk$H, min(mk$Hi), min(hij), gt$cr, gt$mmr, gt$cs, gt$perfect),
    reading = c(psw_h_label(mk$H),
                if (min(mk$Hi) < 0.3) "below .30: this item does not fit the scale" else ".30 or more for every item",
                if (min(hij) < 0) "negative: at least one pair of items orders people in opposite ways" else "no negative pairs",
                if (gt$cr >= 0.9) ".90 or more (conventional criterion)" else "below .90",
                "reproducibility expected from the item marginals alone",
                if (gt$cs >= 0.6) ".60 or more (conventional criterion)" else "below .60",
                sprintf("%d item steps, %s", gt$nsteps,
                        if (gt$order == "entered") "tested in the order entered" else "sorted from easiest to hardest")),
    stringsAsFactors = FALSE)
  if (order == "entered") {
    lo <- apply(X, 2L, min)
    endorsed <- colMeans(sweep(sweep(X, 2L, lo), 2L, pmax(apply(X, 2L, max) - lo, 1e-12), `/`))
    tau <- suppressWarnings(stats::cor(seq_len(ncol(X)), -endorsed, method = "kendall"))
    tab <- rbind(tab, data.frame(
      statistic = "Entered order against observed difficulty (Kendall's tau)", value = tau,
      reading = if (isTRUE(tau > 0.999)) "the items were entered from easiest to hardest"
                else "below 1: some items are easier or harder than their position says",
      stringsAsFactors = FALSE))
  }
  list(table = tab, scalogram = psw_scalogram(X, order))
}

# ---- models ------------------------------------------------------------------

psw_scale_steps <- function(type) {
  if (type == "ordinal") c("congeneric", "tau") else c("congeneric", "tau", "parallel")
}

psw_scale_labels <- function(steps) {
  unname(c(congeneric = "Congeneric", tau = "Tau-equivalent", parallel = "Parallel")[steps])
}

# Single-factor syntax: congeneric (free loadings), tau-equivalent (one common
# loading), parallel (one common loading and one common residual variance).
psw_scale_syntax <- function(prep, step) {
  ids <- prep$items$id
  load <- if (step == "congeneric") ids else paste0("L*", ids)
  out <- paste0("f1 =~ ", paste(load, collapse = " + "))
  if (step == "parallel") out <- c(out, paste0(ids, " ~~ E*", ids))
  paste(out, collapse = "\n")
}

psw_h_label <- function(h) {
  if (is.na(h)) return(NA_character_)
  if (h < 0.3) "not scalable (< .30)" else if (h < 0.4) "weak (.30 to .40)"
  else if (h < 0.5) "medium (.40 to .50)" else "strong (.50 or more)"
}

#' Check one scale: reliability, tau-equivalence and scalability
#'
#' For a set of items meant to measure one construct:
#'
#' * **Reliability**: Cronbach's alpha, Guttman's lambda-2 and omega for the
#'   sum score, plus AVE and the mean inter-item correlation. With ordinal
#'   items, ordinal alpha and ordinal omega (which describe the latent
#'   responses, not the sum score) are added.
#' * **Measurement model**: congeneric, tau-equivalent (equal loadings) and
#'   parallel (equal loadings and residual variances) one-factor models,
#'   each compared with the previous one. Alpha equals reliability only
#'   under tau-equivalence, and unit-weighted sum scores are optimal only
#'   for parallel items, so this is the check behind two routine habits.
#' * **Items**: loading, R-squared, corrected item-rest correlation, alpha
#'   without the item, and the item scalability H.
#' * **Scalability**: Loevinger's H (Mokken scale analysis; works for binary
#'   and for ordinal items and is the probabilistic counterpart of a Guttman
#'   scale) and classical Guttman scalogram statistics computed on item
#'   steps.
#'
#' Classical statistics use the cases with complete data on all items.
#'
#' @inheritParams invarianceWiz
#' @param items Character vector of item names (at least three).
#' @param name Name of the scale, used in tables and plots.
#' @param guttman.order Order of the items the Guttman statistics are tested
#'   against: `"entered"` (default: the order of `items`, easiest first, i.e.
#'   your hypothesis; binary items only, steps of polytomous items are always
#'   sorted) or `"difficulty"` (sorted from the most to the least endorsed).
#'   Loevinger's H does not depend on the order.
#' @param estimator `"auto"` (`"MLR"` for continuous items, `"WLSMV"` for
#'   ordinal items) or one accepted by [invarianceWiz()].
#' @return An object of class `psw_scale`: a list with data frames
#'   `reliability`, `fit`, `comparison` (with verdicts), `items`,
#'   `scalability`, `loadings`, plus `syntax`, `notes` and `info`.
#' @references
#' Graham, J. M. (2006). Congeneric and (essentially) tau-equivalent
#' estimates of score reliability. *Educational and Psychological
#' Measurement, 66*(6), 930-944.
#'
#' Mokken, R. J. (1971). *A theory and procedure of scale analysis*. Mouton.
#'
#' Molenaar, I. W. (1991). A weighted Loevinger H-coefficient extending
#' Mokken scaling to multicategory items. *Kwantitatieve Methoden, 37*,
#' 97-117.
#' @examples
#' \donttest{
#' data(wizdata)
#' scaleWiz(wizdata, paste0("trust", 1:4), name = "Trust")
#' }
#' @export
scaleWiz <- function(data, items, name = "Scale",
                     type = c("auto", "ordinal", "continuous"), estimator = "auto",
                     missing = c("auto", "listwise", "fiml", "pairwise"),
                     criteria = c("cfi", "rmsea"), cutoffs = list(),
                     guttman.order = c("entered", "difficulty"),
                     parallel = "auto", ncores = NULL, keep.fits = FALSE) {
  type <- match.arg(type)
  missing <- match.arg(missing)
  guttman.order <- match.arg(guttman.order)
  factors <- stats::setNames(list(unique(as.character(items))), name)
  prep <- psw_prepare(data, factors, type = type)
  ord <- prep$type == "ordinal"
  notes <- prep$notes
  num <- if (ord) psw_prepare(data, factors, type = "continuous") else prep

  estimator <- toupper(estimator)
  if (estimator == "AUTO") estimator <- if (ord) "WLSMV" else "MLR"
  allowed <- if (ord) c("WLSMV", "ULSMV", "DWLS", "ULS") else psw_ml_estimators
  if (!estimator %in% allowed)
    stop(sprintf("Estimator '%s' is not available for %s items. Choose one of: %s.",
                 estimator, prep$type, paste(allowed, collapse = ", ")), call. = FALSE)
  miss <- psw_resolve_missing(missing, estimator, ord, prep$has_missing)
  if (ord && miss == "fiml") stop("FIML is not available for ordinal items.", call. = FALSE)
  if (!ord && miss == "pairwise") stop("Pairwise deletion is only available for ordinal items.", call. = FALSE)

  steps <- psw_scale_steps(prep$type)
  labels <- psw_scale_labels(steps)
  syntax <- stats::setNames(lapply(steps, function(s) psw_scale_syntax(prep, s)), steps)
  tasks <- lapply(steps, function(s)
    list(kind = "inv", step = s, model = syntax[[s]], prep = prep, estimator = estimator,
         missing = miss, parameterization = "delta", culprits = FALSE))
  # with ordinal items, omega for the sum score needs the items-as-scored model too
  if (ord)
    tasks <- c(tasks, list(list(kind = "inv", step = "scored", model = psw_scale_syntax(num, "congeneric"),
                                prep = num, estimator = "MLR",
                                missing = if (num$has_missing) "fiml" else "listwise",
                                parameterization = "delta", culprits = FALSE)))

  t0 <- proc.time()[["elapsed"]]
  run <- psw_run_tasks(tasks, parallel = parallel, ncores = ncores)
  res <- stats::setNames(run$results[seq_along(steps)], steps)
  scored <- if (ord) run$results[[length(tasks)]] else res$congeneric
  lad <- psw_ladder_tables(res, steps, labels, ord)
  notes <- c(notes, lad$notes)

  # ---- reliability ----
  X <- as.matrix(num$data[num$items$id])
  cc <- stats::complete.cases(X)
  X <- X[cc, , drop = FALSE]
  colnames(X) <- num$items$name
  k <- ncol(X)
  if (nrow(X) < 10L) stop("Fewer than 10 cases have complete data on all items.", call. = FALSE)
  C <- stats::cov(X)
  R <- stats::cov2cor(C)
  rel <- list()
  add <- function(coefficient, value, basis)
    rel[[length(rel) + 1L]] <<- data.frame(coefficient = coefficient, value = value, basis = basis,
                                           stringsAsFactors = FALSE)
  add("Cronbach's alpha", psw_alpha(C), "sum score")
  add("Guttman's lambda-2", psw_lambda2(C), "sum score")
  omega_of <- function(r) {
    if (!isTRUE(r$converged)) return(c(NA_real_, NA_real_))
    est <- lavaan::lavInspect(r$fit, "est")
    lam <- as.numeric(est$lambda)
    th <- diag(est$theta)
    std <- r$loadings$std
    c(sum(lam)^2 / (sum(lam)^2 + sum(th)), mean(std^2))
  }
  om <- omega_of(scored)
  add("McDonald's omega", om[1L], "sum score")
  if (ord && isTRUE(res$congeneric$converged)) {
    P <- unclass(lavaan::lavInspect(res$congeneric$fit, "sampstat")$cov)
    add("Ordinal alpha", psw_alpha(P), "latent responses")
    l <- res$congeneric$loadings$std
    add("Ordinal omega", sum(l)^2 / (sum(l)^2 + sum(1 - l^2)), "latent responses")
    add("AVE", mean(l^2), "latent responses")
  } else {
    add("AVE", om[2L], "items")
  }
  add("Mean inter-item correlation", mean(R[lower.tri(R)]), "items")
  reliability <- do.call(rbind, rel)

  # ---- items ----
  mk <- psw_mokken_h(X)
  load <- res$congeneric$loadings
  items_df <- data.frame(
    item = num$items$name,
    loading = if (is.null(load)) NA_real_ else load$std[match(num$items$name, load$item)],
    mean = colMeans(X), sd = apply(X, 2L, stats::sd),
    rest.r = vapply(seq_len(k), function(j) stats::cor(X[, j], rowSums(X[, -j, drop = FALSE])), numeric(1)),
    alpha.drop = vapply(seq_len(k), function(j) psw_alpha(C[-j, -j, drop = FALSE]), numeric(1)),
    H = unname(mk$Hi), stringsAsFactors = FALSE, row.names = NULL)
  items_df$r2 <- items_df$loading^2

  # ---- scalability, for both item orders (cheap), so jamovi can switch without refitting ----
  guttman <- list(entered = psw_scalability(X, mk, "entered"), difficulty = psw_scalability(X, mk, "difficulty"))
  scalability <- guttman[[guttman.order]]$table
  if (any(apply(X, 2L, function(v) length(unique(v))) > 2L))
    notes <- c(notes, paste0("Guttman statistics are computed on item steps (response at or above each category), sorted from easiest to hardest. ",
                             "The steps of one item are nested by construction, which inflates CR; ",
                             "read CS and Loevinger's H instead."))
  if (sum(!cc) > 0L)
    notes <- c(notes, sprintf("Reliability, item and scalability statistics use the %d cases with complete data (%d dropped).",
                              nrow(X), sum(!cc)))

  out <- list(
    reliability = reliability, fit = lad$fit, comparison = lad$comparison, items = items_df,
    scalability = scalability, guttman = guttman, hij = mk$Hij,
    loadings = if (is.null(load)) NULL else cbind(step = "congeneric", label = "Congeneric", load, stringsAsFactors = FALSE),
    cors = NULL,
    syntax = lapply(syntax, psw_syntax_display, prep = prep), notes = unique(notes),
    info = list(name = name, type = prep$type, estimator = estimator, missing = miss,
                variant = lad$variant, factors = factors, n = nrow(prep$data), n.complete = nrow(X),
                mode = run$mode, workers = run$workers, time = proc.time()[["elapsed"]] - t0,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- lapply(res, function(r) r$fit)
  class(out) <- "psw_scale"
  invarianceVerdict(out, criteria = criteria, cutoffs = cutoffs)
}

# What the model comparison means for everyday practice, in one sentence each.
psw_scale_reading <- function(x) {
  cp <- x$comparison
  v <- stats::setNames(cp$verdict, cp$step)
  out <- character()
  if ("tau" %in% names(v))
    out <- c(out, switch(v[["tau"]],
      "Supported" = "Equal loadings are tenable: alpha is an adequate estimate of reliability for this scale.",
      "Not supported" = "Loadings differ: alpha underestimates reliability here; report omega instead.",
      "Ambiguous" = "The evidence on equal loadings is mixed: report omega alongside alpha.",
      "Equal loadings could not be evaluated."))
  if ("parallel" %in% names(v))
    out <- c(out, switch(v[["parallel"]],
      "Supported" = "Equal residual variances are tenable too: the items are interchangeable and unit weights are optimal.",
      "Not supported" = "Residual variances differ: items are not interchangeable, and weighted scores would be slightly more reliable than the plain sum.",
      "Ambiguous" = "The evidence on equal residual variances is mixed.",
      "Parallelism could not be evaluated."))
  out
}

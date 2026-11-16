# engine-cfa.R -- one plain CFA (one factor, or several correlated factors) in
# one group: fit, loadings, factor correlations, reliability, the largest
# modification indices and the factor correlations under other ways of scoring.
#
# This is the entry-level analysis. The others start from the same model and
# ask one more question each: across groups (engine-invariance.R), across
# estimators (engine-multiest.R), with cross-loadings (engine-esem.R).

# Estimator and missing-data handling for one single-group model, with the
# checks every analysis needs. Returns list(estimator, missing).
psw_resolve_estimation <- function(estimator, missing, prep) {
  ord <- prep$type == "ordinal"
  estimator <- toupper(estimator)
  if (estimator == "AUTO") estimator <- if (ord) "WLSMV" else "MLR"
  allowed <- if (ord) c("WLSMV", "ULSMV", "DWLS", "ULS") else psw_ml_estimators
  if (!estimator %in% allowed)
    stop(sprintf("Estimator '%s' is not available for %s items. Choose one of: %s.",
                 estimator, prep$type, paste(allowed, collapse = ", ")), call. = FALSE)
  miss <- psw_resolve_missing(missing, estimator, ord, prep$has_missing)
  if (ord && miss == "fiml") stop("FIML is not available for ordinal items; use 'pairwise' or 'listwise'.", call. = FALSE)
  if (!ord && miss == "pairwise") stop("Pairwise deletion is only available for ordinal items.", call. = FALSE)
  if (miss == "fiml" && estimator == "MLM") stop("FIML is not available with MLM; use ML or MLR.", call. = FALSE)
  list(estimator = estimator, missing = miss)
}

# Correlations between factor scores estimated from a fitted model: regression
# scores, and Bartlett scores (continuous items only). NULL when unavailable.
psw_score_cors <- function(fit, ordinal) {
  one <- function(...) tryCatch(
    suppressWarnings(stats::cor(lavaan::lavPredict(fit, ...), use = "pairwise.complete.obs")),
    error = function(e) NULL)
  list(scores = one(), bartlett = if (ordinal) NULL else one(method = "Bartlett"))
}

# Composite (item mean) scores per factor: their correlations, each scale's
# alpha, and the correlations corrected for unreliability with those alphas.
# Items are always used as scored (numeric), also when the model treats them
# as ordinal, because that is what a researcher computing mean scores does.
psw_composite_cors <- function(data, factors) {
  num <- psw_prepare(data, factors, type = "continuous")
  fids <- names(num$fnames)
  ids_of <- function(f) num$items$id[num$items$fid == f]
  comp <- vapply(fids, function(f) rowMeans(num$data[ids_of(f)], na.rm = TRUE), numeric(nrow(num$data)))
  alpha <- vapply(fids, function(f) psw_alpha(stats::cov(num$data[ids_of(f)], use = "complete.obs")), numeric(1))
  r <- stats::cor(comp, use = "pairwise.complete.obs")
  list(composite = r, corrected = r / sqrt(outer(alpha, alpha)), alpha = unname(alpha))
}

# The largest modification indices of a single-group CFA, of the two kinds a
# reader would act on: a residual correlation between two items, or an item
# loading on another factor. NULL if lavaan cannot compute them.
psw_suggestions <- function(fit, prep) {
  it <- prep$items
  mi <- tryCatch(suppressWarnings(lavaan::modindices(fit, sort. = FALSE)), error = function(e) NULL)
  if (is.null(mi) || !nrow(mi)) return(NULL)
  res <- mi$op == "~~" & mi$lhs %in% it$id & mi$rhs %in% it$id & mi$lhs != mi$rhs
  cross <- mi$op == "=~" & mi$lhs %in% names(prep$fnames) & mi$rhs %in% it$id
  mi <- mi[(res | cross) & is.finite(mi$mi), , drop = FALSE]
  if (!nrow(mi)) return(NULL)
  is_res <- mi$op == "~~"
  name <- function(id) ifelse(id %in% it$id, it$name[match(id, it$id)], unname(prep$fnames[id]))
  same <- it$fid[match(mi$lhs, it$id)] == it$fid[match(mi$rhs, it$id)]
  out <- data.frame(
    kind = ifelse(is_res, "Residual correlation", "Cross-loading"),
    lhs = name(mi$lhs), rhs = name(mi$rhs),
    detail = ifelse(is_res, ifelse(same, "same factor", "different factors"), ""),
    mi = mi$mi, epc = mi$sepc.all, stringsAsFactors = FALSE)
  out <- out[order(-out$mi), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Confirmatory factor analysis: one factor or several correlated factors
#'
#' Fits one CFA in one group and returns what a report usually needs: fit
#' indices, standardised loadings, factor correlations and reliability per
#' factor. Items are treated as ordinal (WLSMV, polychoric correlations) or as
#' continuous (MLR) according to their type.
#'
#' Two extras come with it:
#'
#' * `suggestions`: the largest modification indices, restricted to residual
#'   correlations between items and cross-loadings. They show where the model
#'   strains; they are not a to-do list.
#' * `correlations` (two or more factors): each factor correlation as it comes
#'   out of composite (item mean) scores, composites corrected for
#'   unreliability, the CFA itself, and CFA factor scores. This is the answer
#'   to "would I have found the same with plain mean scores?".
#'
#' Reliability: alpha is computed from the items as scored. Omega is the
#' composite reliability implied by the model; for continuous items it refers
#' to the sum score, for ordinal items to the latent responses (so it is
#' usually higher than what the sum score achieves).
#'
#' Two additions to the model are possible. `residuals` lets the unique parts
#' of two items correlate (for example two items with nearly the same
#' wording). `higher` replaces the correlations among three or more factors
#' by one second-order factor on which they all load.
#'
#' @inheritParams invarianceWiz
#' @param residuals Residual correlations to add: a list of pairs of item
#'   names, e.g. `list(c("x1", "x2"))`.
#' @param higher Name of a second-order factor, or `NULL` (default) for
#'   correlated factors. Needs at least three factors.
#' @return An object of class `psw_cfa`: a list with data frames `fit` (one
#'   row), `loadings`, `cors`, `reliability`, `suggestions`, `correlations`,
#'   `residuals`, the list `higher` (name and loadings of the factors on the
#'   second-order factor), plus `syntax`, `notes` and `info`.
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' factorWiz(wizdata, f)
#' }
#' @export
factorWiz <- function(data, factors, type = c("auto", "ordinal", "continuous"), estimator = "auto",
                      missing = c("auto", "listwise", "fiml", "pairwise"),
                      residuals = NULL, higher = NULL, keep.fits = FALSE) {
  type <- match.arg(type)
  missing <- match.arg(missing)
  prep <- psw_prepare(data, factors, type = type)
  ord <- prep$type == "ordinal"
  est <- psw_resolve_estimation(estimator, missing, prep)
  it <- prep$items
  fn <- unname(prep$fnames)
  notes <- prep$notes

  fids <- names(prep$fnames)
  syntax <- psw_syntax(prep)
  # residual correlations: pairs of item names -> "i3 ~~ i5"
  pairs <- unique(lapply(residuals, function(p) sort(as.character(unlist(p)))))
  pairs <- Filter(function(p) length(p) == 2L && p[1L] != p[2L], pairs)
  bad <- setdiff(unlist(pairs), it$name)
  if (length(bad)) stop("Residual correlations: not an item of the model: ", paste(bad, collapse = ", "), call. = FALSE)
  if (length(pairs))
    syntax <- paste(c(syntax, vapply(pairs, function(p) paste(it$id[match(p, it$name)], collapse = " ~~ "), "")), collapse = "\n")
  # second-order factor "h" instead of factor correlations
  if (!is.null(higher)) {
    if (length(fids) < 3L)
      stop("A second-order factor needs at least three factors (with two it is the same model as two correlated factors).", call. = FALSE)
    syntax <- paste(c(syntax, paste0("h =~ ", paste(fids, collapse = " + "))), collapse = "\n")
  }
  r <- psw_fit(syntax, prep, est$estimator, est$missing)
  if (!is.null(r$error)) stop(r$error, call. = FALSE)
  if (!psw_converged(r$fit)) stop("The model did not converge.", call. = FALSE)
  notes <- c(notes, r$warnings)

  fm <- psw_fitmeasures_raw(r$fit)
  variant <- psw_choose_variant(list(fm), ord)
  fit <- data.frame(label = "CFA", t(psw_measures(fm, variant)), stringsAsFactors = FALSE)
  if (isTRUE(fit$df == 0))
    notes <- c(notes, "The model has zero degrees of freedom (one factor with three items): it fits perfectly by construction, so the fit indices say nothing.")

  e <- psw_estimates(r$fit, prep)
  loadings <- e$loadings
  loadings$r2 <- loadings$std^2

  pe <- lavaan::parameterEstimates(r$fit, standardized = TRUE, ci = FALSE)
  rs <- pe[pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% it$id & pe$rhs %in% it$id, , drop = FALSE]
  resid <- if (nrow(rs)) data.frame(item1 = it$name[match(rs$lhs, it$id)], item2 = it$name[match(rs$rhs, it$id)],
                                    r = rs$std.all, p = rs$pvalue, stringsAsFactors = FALSE)
  hi <- pe[pe$op == "=~" & pe$lhs == "h", , drop = FALSE]
  higher_out <- if (!is.null(higher))
    list(name = higher, loadings = data.frame(factor = unname(prep$fnames[hi$rhs]), est = hi$est, se = hi$se,
                                              p = hi$pvalue, std = hi$std.all, stringsAsFactors = FALSE))

  # ---- reliability per factor ----
  cc <- psw_composite_cors(data, factors)
  theta <- unclass(lavaan::lavInspect(r$fit, "est")$theta)
  phi <- diag(unclass(lavaan::lavInspect(r$fit, "cov.lv")))
  reliability <- do.call(rbind, lapply(seq_along(fn), function(f) {
    ids <- it$id[it$fid == fids[f]]
    l <- loadings[loadings$factor == fn[f], , drop = FALSE]
    # omega of the sum score from the raw solution: true-score variance is
    # (sum of loadings)^2 x factor variance, error variance is the sum of the
    # residual (co)variances of the items. For ordinal items the raw solution
    # is on the latent-response scale.
    true <- sum(l$est)^2 * phi[[fids[f]]]
    data.frame(factor = fn[f], k = length(ids), alpha = cc$alpha[f],
               omega = true / (true + sum(theta[ids, ids])), ave = mean(l$std^2), stringsAsFactors = FALSE)
  }))

  # ---- factor correlations by scoring method ----
  correlations <- NULL
  if (length(fn) >= 2L) {
    pairs <- which(upper.tri(diag(length(fn))), arr.ind = TRUE)
    pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
    pick <- function(M) if (is.null(M)) rep(NA_real_, nrow(pairs)) else unname(M[pairs])
    sc <- lapply(psw_score_cors(r$fit, ord), function(M) if (!is.null(M)) M[fids, fids])
    correlations <- data.frame(
      factor1 = fn[pairs[, 1L]], factor2 = fn[pairs[, 2L]],
      composite = pick(cc$composite), corrected = pick(cc$corrected),
      cfa = pick(unclass(lavaan::lavInspect(r$fit, "cor.lv"))[fids, fids]),
      scores = pick(sc$scores), bartlett = pick(sc$bartlett), stringsAsFactors = FALSE)
    if (any(cc$alpha <= 0, na.rm = TRUE))
      notes <- c(notes, "A scale has alpha at or below zero, so the corrected composite correlation is not meaningful.")
  }

  out <- list(
    fit = fit, loadings = loadings, cors = e$cors, reliability = reliability,
    suggestions = psw_suggestions(r$fit, prep), correlations = correlations,
    residuals = resid, higher = higher_out,
    syntax = psw_syntax_display(syntax, prep), notes = unique(notes),
    info = list(type = prep$type, estimator = est$estimator, missing = est$missing, variant = variant,
                factors = split(it$name, factor(it$factor, levels = fn)),
                n = nrow(prep$data), mode = "sequential", workers = 0L, time = r$time,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- list(cfa = r$fit)
  class(out) <- "psw_cfa"
  out
}

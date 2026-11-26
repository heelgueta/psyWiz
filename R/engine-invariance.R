# engine-invariance.R -- the measurement invariance ladder.

# One step of the ladder; runs in the main process or in a worker.
psw_inv_task <- function(task) {
  r <- psw_fit(task$model, task$prep, task$estimator, task$missing, task$parameterization)
  r$step <- task$step
  r$converged <- psw_converged(r$fit)
  if (is.null(r$error) && !r$converged) r$error <- "The model did not converge."
  r$fm <- psw_fitmeasures_raw(r$fit)
  est <- psw_estimates(r$fit, task$prep)
  r$loadings <- est$loadings
  r$cors <- est$cors
  r$means <- est$means
  r$culprits <- if (r$converged && task$step != "configural" && !isFALSE(task$culprits))
    psw_culprits(r$fit, task$step, task$prep) else NULL
  r
}

#' Which constraints of a step fit worst? (internal)
#'
#' Univariate score (Lagrange multiplier) tests for the constraints a step
#' adds: how much the model chi-square is expected to drop if that single
#' constraint were released. For equality constraints this is
#' `lavaan::lavTestScore()`; for parameters fixed to a constant in the
#' non-reference groups (ordinal intercepts and residual variances) it is the
#' modification index.
#' @return Data frame (item, factor, detail, group, x2, df, p) sorted by x2, or NULL.
#' @keywords internal
#' @noRd
psw_culprits <- function(fit, step, prep) {
  it <- prep$items
  ord <- prep$type == "ordinal"
  res <- tryCatch({
    if (ord && step %in% c("scalar", "strict")) {
      op <- if (step == "scalar") "~1" else "~~"
      mi <- suppressWarnings(lavaan::modindices(fit, op = op, sort. = FALSE))
      if (is.null(mi$group)) mi$group <- 1L
      keep <- mi$lhs %in% it$id & mi$group > 1L
      if (op == "~~") keep <- keep & mi$lhs == mi$rhs
      mi <- mi[keep, , drop = FALSE]
      data.frame(id = mi$lhs, detail = rep("", nrow(mi)), group = mi$group, x2 = mi$mi,
                 stringsAsFactors = FALSE)
    } else {
      pt <- as.data.frame(lavaan::parTable(fit), stringsAsFactors = FALSE)
      eq <- pt[pt$op == "==", , drop = FALSE]
      a <- pt[match(eq$lhs, pt$plabel), , drop = FALSE]
      b <- pt[match(eq$rhs, pt$plabel), , drop = FALSE]
      want <- c(thresholds = "|", metric = "=~", scalar = "~1", strict = "~~")[[step]]
      keep <- !is.na(b$op) & b$op == want
      if (want == "~~") keep <- keep & b$lhs == b$rhs
      idx <- which(keep)
      x2 <- psw_score_uni(fit, idx)
      b <- b[idx, , drop = FALSE]
      data.frame(id = if (want == "=~") b$rhs else b$lhs,
                 detail = if (want == "|") b$rhs else rep("", nrow(b)),
                 group = pmax(a$group[idx], b$group), x2 = x2,
                 stringsAsFactors = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  res <- res[res$id %in% it$id & is.finite(res$x2), , drop = FALSE]
  if (nrow(res) == 0L) return(NULL)
  m <- match(res$id, it$id)
  out <- data.frame(step = step, parameter = psw_step_target(step),
                    item = it$name[m], factor = it$factor[m], detail = res$detail,
                    group = prep$groups[res$group], x2 = res$x2, df = 1L,
                    p = stats::pchisq(res$x2, 1, lower.tail = FALSE),
                    stringsAsFactors = FALSE)
  out <- out[order(-out$x2), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Univariate score tests for releasing single equality constraints.
#
# Gives the same numbers as lavaan::lavTestScore(fit, release = idx)$uni$X2,
# which inverts a bordered information matrix once per constraint (minutes for
# models with hundreds of constraints). Here the matrix
#   Z = [ I  R' ]      I = expected information, R = constraint Jacobian
#       [ R  0  ]
# is inverted once. With Z^-1 = [A B'; B C], the Lagrange multipliers are
# lambda = B s (s = score) and the test for releasing constraint r alone is
# N * lambda_r^2 / (-C_rr), which follows from the formula for the inverse of
# a matrix with one row and column removed.
psw_score_uni <- function(fit, idx) {
  R <- unclass(fit@Model@con.jac)[, , drop = FALSE]
  score <- as.numeric(lavaan::lavTech(fit, "gradient.logl"))
  info <- unclass(lavaan::lavTech(fit, "information.expected"))
  N <- lavaan::lavInspect(fit, "ntotal")
  np <- nrow(info)
  nc <- nrow(R)
  Z <- rbind(cbind(info, t(R)), cbind(R, matrix(0, nc, nc)))
  Zi <- tryCatch(solve(Z), error = function(e) MASS::ginv(Z))
  con <- np + seq_len(nc)
  lambda <- as.numeric(Zi[con, seq_len(np), drop = FALSE] %*% score)
  v <- -diag(Zi)[con]
  (N * lambda^2 / v)[idx]
}

# Chi-square difference test between two nested fits (fit0 less constrained).
psw_lrt <- function(fit0, fit1) {
  na <- list(dchisq = NA_real_, ddf = NA_real_, p = NA_real_, note = NULL)
  df0 <- lavaan::fitMeasures(fit0, "df")[[1L]]
  df1 <- lavaan::fitMeasures(fit1, "df")[[1L]]
  if (df1 - df0 < 1) {
    na$ddf <- df1 - df0
    return(na)
  }
  warns <- character()
  res <- tryCatch(
    withCallingHandlers(lavaan::lavTestLRT(fit0, fit1),
                        warning = function(w) {
                          warns <<- c(warns, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) e)
  if (inherits(res, "error")) {
    na$ddf <- df1 - df0
    na$note <- psw_clean_msg(conditionMessage(res))
    return(na)
  }
  list(dchisq = res[2L, "Chisq diff"], ddf = res[2L, "Df diff"], p = res[2L, "Pr(>Chisq)"],
       note = if (length(warns)) paste(unique(psw_clean_msg(warns)), collapse = " ") else NULL)
}

#' Fit and step-comparison tables for a sequence of nested models (internal)
#'
#' Shared by the invariance ladder and the congeneric / tau-equivalent /
#' parallel ladder: one row of fit indices per model, and each model compared
#' with the previous one that converged (chi-square difference test and
#' changes in CFI, RMSEA and SRMR).
#' @param res Named list of task results (each with `fit`, `fm`, `converged`).
#' @return list(fit, comparison, variant, notes)
#' @keywords internal
#' @noRd
psw_ladder_tables <- function(res, steps, labels, ordinal) {
  notes <- character()
  variant <- psw_choose_variant(lapply(res, function(r) r$fm), ordinal)
  fit <- do.call(rbind, lapply(seq_along(steps), function(i) {
    r <- res[[i]]
    data.frame(step = steps[i], label = labels[i], converged = r$converged,
               t(psw_measures(r$fm, variant)), time = r$time, stringsAsFactors = FALSE)
  }))
  rownames(fit) <- NULL
  for (i in seq_along(steps)) {
    r <- res[[i]]
    if (!is.null(r$error)) notes <- c(notes, sprintf("%s model: %s", labels[i], r$error))
    if (length(r$warnings)) notes <- c(notes, sprintf("%s model: %s", labels[i], r$warnings))
  }
  comparison <- do.call(rbind, lapply(seq_along(steps)[-1L], function(i) {
    row <- data.frame(step = steps[i], label = labels[i], vs = NA_character_,
                      dchisq = NA_real_, ddf = NA_real_, p = NA_real_,
                      dcfi = NA_real_, drmsea = NA_real_, dsrmr = NA_real_,
                      stringsAsFactors = FALSE)
    prev <- rev(which(fit$converged[seq_len(i - 1L)]))[1L]
    if (!fit$converged[i] || is.na(prev)) return(row)
    row$vs <- labels[prev]
    lrt <- psw_lrt(res[[prev]]$fit, res[[i]]$fit)
    row$dchisq <- lrt$dchisq
    row$ddf <- lrt$ddf
    row$p <- lrt$p
    if (!is.null(lrt$note)) notes <<- c(notes, sprintf("%s vs %s: %s", labels[i], labels[prev], lrt$note))
    row$dcfi <- fit$cfi[i] - fit$cfi[prev]
    row$drmsea <- fit$rmsea[i] - fit$rmsea[prev]
    row$dsrmr <- fit$srmr[i] - fit$srmr[prev]
    row
  }))
  rownames(comparison) <- NULL
  list(fit = fit, comparison = comparison, variant = variant, notes = notes)
}

#' Measurement invariance across groups
#'
#' Fits a one-factor or correlated-factors CFA in several groups under
#' increasingly strict equality constraints and compares each step with the
#' previous one.
#'
#' For continuous items the steps are configural, metric (loadings), scalar
#' (loadings and intercepts) and, optionally, strict (residual variances). For
#' ordinal items the ladder of Wu and Estabrook (2016) is used: configural,
#' thresholds, thresholds and loadings, thresholds, loadings and intercepts.
#' The models and decision rules are described in `docs/METHODS.md` in the
#' package repository.
#'
#' @param data A data frame.
#' @param factors A named list giving the items of each factor, for example
#'   `list(Verbal = c("v1", "v2", "v3"), Spatial = c("s1", "s2", "s3"))`. A
#'   plain character vector is read as a single factor.
#' @param group Name of the grouping variable.
#' @param type How to treat the items: `"auto"` (ordinal if they are factors,
#'   continuous if numeric), `"ordinal"` or `"continuous"`.
#' @param estimator `"auto"` (`"MLR"` for continuous items, `"WLSMV"` for
#'   ordinal items), or one of `"ML"`, `"MLR"`, `"MLM"` (continuous) or
#'   `"WLSMV"`, `"ULSMV"`, `"DWLS"`, `"ULS"` (ordinal).
#' @param missing `"auto"`, `"listwise"`, `"fiml"` or `"pairwise"`. `"auto"`
#'   uses full-information ML for ML/MLR and pairwise deletion for ordinal
#'   items when there are missing values.
#' @param strict Also test equality of residual variances?
#' @param parameterization `"auto"`, `"delta"` or `"theta"` (ordinal items
#'   only). The strict step needs `"theta"`, which `"auto"` selects for you.
#' @param criteria Which criteria decide the verdict; any of `"cfi"`,
#'   `"rmsea"`, `"srmr"`, `"chisq"`.
#' @param cutoffs Named list overriding the default cutoffs, see
#'   [invarianceVerdict()].
#' @param parallel `"auto"`, `TRUE` or `FALSE`: fit the steps on several
#'   cores? `"auto"` does so only when the first model is slow.
#' @param ncores Maximum number of worker processes.
#' @param keep.fits Keep the fitted lavaan objects in the result (`$fits`)?
#' @param partial Also run [invariancePartial()] and store its result in
#'   `$partial`?
#' @param max.free,min.invariant Passed to [invariancePartial()].
#'
#' @return An object of class `psw_invariance`: a list with data frames `fit`
#'   (one row per step), `comparison` (each step against the previous one,
#'   with a verdict), `culprits` (score tests for the constraints of each
#'   step), `loadings`, `groups`, plus `syntax`, `notes`, `info`, and optionally
#'   `partial` (if `partial = TRUE`).
#'
#' @references
#' Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
#' measurement invariance. *Structural Equation Modeling, 14*(3), 464-504.
#'
#' Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
#' analysis models of different levels of invariance for ordered categorical
#' outcomes. *Psychometrika, 81*(4), 1014-1045.
#'
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' res <- invarianceWiz(wizdata, f, group = "country")
#' res
#' }
#' @export
invarianceWiz <- function(data, factors, group,
                          type = c("auto", "ordinal", "continuous"),
                          estimator = "auto",
                          missing = c("auto", "listwise", "fiml", "pairwise"),
                          strict = FALSE,
                          parameterization = c("auto", "delta", "theta"),
                          criteria = c("cfi", "rmsea"), cutoffs = list(),
                          parallel = "auto", ncores = NULL, keep.fits = FALSE,
                          partial = FALSE, max.free = 2L, min.invariant = 2L,
                          means = FALSE, means.ref = NULL) {
  type <- match.arg(type)
  missing <- match.arg(missing)
  parameterization <- match.arg(parameterization)
  if (base::missing(group) || is.null(group) || length(group) != 1L)
    stop("`group` must name one grouping variable.", call. = FALSE)

  prep <- psw_prepare(data, factors, group = group, type = type)
  ord <- prep$type == "ordinal"
  notes <- prep$notes

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

  if (!ord) parameterization <- NA_character_
  else if (parameterization == "auto") parameterization <- if (strict) "theta" else "delta"
  else if (strict && parameterization == "delta")
    stop("The strict step for ordinal items needs parameterization = 'theta'.", call. = FALSE)
  par_arg <- if (ord) parameterization else "delta"

  if (ord && any(prep$nthr == 1L))
    notes <- c(notes, paste0("Binary items: their residual variances are not identified, so differences in ",
                             "thresholds, loadings and residual variances cannot be fully separated ",
                             "(Wu & Estabrook, 2016)."))

  steps <- psw_steps(prep$type, strict)
  labels <- psw_step_labels(steps, prep$type)
  syntax <- stats::setNames(lapply(steps, function(s) psw_syntax(prep, s, par_arg)), steps)
  tasks <- lapply(steps, function(s)
    list(kind = "inv", step = s, model = syntax[[s]], prep = prep, estimator = estimator,
         missing = miss, parameterization = par_arg))

  t0 <- proc.time()[["elapsed"]]
  run <- psw_run_tasks(tasks, parallel = parallel, ncores = ncores)
  res <- stats::setNames(run$results, steps)

  lad <- psw_ladder_tables(res, steps, labels, ord)
  fit <- lad$fit
  comparison <- lad$comparison
  variant <- lad$variant
  notes <- c(notes, lad$notes)

  culprits <- do.call(rbind, lapply(res, function(r) r$culprits))
  if (!is.null(culprits)) {
    culprits$label <- labels[match(culprits$step, steps)]
    rownames(culprits) <- NULL
  }
  loadings <- do.call(rbind, lapply(seq_along(steps), function(i) {
    l <- res[[i]]$loadings
    if (is.null(l) || nrow(l) == 0L) return(NULL)
    cbind(step = steps[i], label = labels[i], l, stringsAsFactors = FALSE)
  }))
  cors <- do.call(rbind, lapply(seq_along(steps), function(i) {
    l <- res[[i]]$cors
    if (is.null(l) || nrow(l) == 0L) return(NULL)
    cbind(step = steps[i], label = labels[i], l, stringsAsFactors = FALSE)
  }))

  out <- list(
    fit = fit, comparison = comparison, culprits = culprits, loadings = loadings, cors = cors,
    groups = data.frame(group = names(prep$n), n = as.integer(prep$n), stringsAsFactors = FALSE),
    syntax = lapply(syntax, psw_syntax_display, prep = prep),
    notes = unique(notes),
    # the prepared data, so that invariancePartial() can refit; the jamovi layer drops it before caching
    prep = prep,
    means_scalar = if (!is.null(res$scalar)) res$scalar$means else NULL,
    means_strict = if (!is.null(res$strict)) res$strict$means else NULL,
    info = list(type = prep$type, estimator = estimator, missing = miss,
                parameterization = parameterization, variant = variant,
                factors = split(prep$items$name, factor(prep$items$factor, levels = unname(prep$fnames))),
                group = group, nthr = stats::setNames(prep$nthr, prep$items$name),
                mode = run$mode, workers = run$workers,
                time = proc.time()[["elapsed"]] - t0,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- lapply(res, function(r) r$fit)
  class(out) <- "psw_invariance"
  out <- invarianceVerdict(out, criteria = criteria, cutoffs = cutoffs)

  if (isTRUE(partial)) {
    out$partial <- invariancePartial(out, auto = TRUE, max.free = max.free,
                                     min.invariant = min.invariant,
                                     criteria = criteria, cutoffs = cutoffs)
  }
  if (isTRUE(means)) {
    out$means <- invarianceMeans(out, reference = means.ref)
  }
  out
}

psw_default_cutoffs <- list(cfi = 0.010, rmsea = 0.015, srmr.loadings = 0.030, srmr = 0.010, alpha = 0.05)

#' Apply decision criteria to an invariance ladder
#'
#' Adds (or recomputes) the flags and the verdict of each step comparison
#' without refitting any model, so cutoffs can be changed cheaply.
#'
#' A step is flagged on a criterion when, relative to the previous step,
#' CFI drops by more than `cutoffs$cfi`, RMSEA rises by more than
#' `cutoffs$rmsea`, SRMR rises by more than `cutoffs$srmr.loadings` (the step
#' that constrains loadings) or `cutoffs$srmr` (other steps), or the
#' chi-square difference test has p < `cutoffs$alpha`. Defaults follow Chen
#' (2007). The verdict is `"Supported"` when none of the chosen criteria is
#' flagged, `"Not supported"` when all of them are, and `"Ambiguous"` when they
#' disagree.
#'
#' @param x A `psw_invariance` object (or the result of [scaleWiz()]).
#' @param criteria Any of `"cfi"`, `"rmsea"`, `"srmr"`, `"chisq"`.
#' @param cutoffs Named list with any of `cfi` (0.010), `rmsea` (0.015),
#'   `srmr.loadings` (0.030), `srmr` (0.010), `alpha` (0.05).
#' @return `x` with columns `flag.cfi`, `flag.rmsea`, `flag.srmr`,
#'   `flag.chisq` and `verdict` in `x$comparison`.
#' @references
#' Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
#' measurement invariance. *Structural Equation Modeling, 14*(3), 464-504.
#' @export
invarianceVerdict <- function(x, criteria = c("cfi", "rmsea"), cutoffs = list()) {
  stopifnot(inherits(x, c("psw_invariance", "psw_scale")))
  if (length(criteria))
    criteria <- match.arg(criteria, c("cfi", "rmsea", "srmr", "chisq"), several.ok = TRUE)
  cut <- utils::modifyList(psw_default_cutoffs, as.list(cutoffs))
  cp <- x$comparison
  eps <- 1e-9
  cp$flag.cfi <- cp$dcfi < -cut$cfi - eps
  cp$flag.rmsea <- cp$drmsea > cut$rmsea + eps
  cp$flag.srmr <- cp$dsrmr > ifelse(cp$step %in% c("metric", "tau"), cut$srmr.loadings, cut$srmr) + eps
  cp$flag.chisq <- cp$p < cut$alpha
  cp$verdict <- vapply(seq_len(nrow(cp)), function(i) {
    if (is.na(cp$vs[i])) return("Not estimated")
    if (!is.na(cp$ddf[i]) && cp$ddf[i] < 1) return("Equivalent to previous model")
    if (!length(criteria)) return("No criterion selected")
    fl <- unlist(cp[i, paste0("flag.", criteria), drop = TRUE])
    fl <- fl[!is.na(fl)]
    if (!length(fl)) return("Not available")
    if (!any(fl)) "Supported" else if (all(fl)) "Not supported" else "Ambiguous"
  }, "")
  x$comparison <- cp
  x$criteria <- criteria
  x$cutoffs <- cut
  x
}

#' Partial measurement invariance
#'
#' Refits the invariance ladder, releasing equality constraints item by item
#' when a step is not supported. It answers the question that follows a
#' failed step: is the problem confined to a few items, and how far does
#' invariance go once those are set free?
#'
#' Search (when `auto = TRUE`). At a step whose verdict is "Not supported" or
#' "Ambiguous", the item with the largest score test among the constraints
#' that step adds is released in every group, the step is refitted and
#' judged again against the same previous model, and so on until the step is
#' supported or no further release is allowed. A release is allowed only if
#' it leaves at least `min.invariant` items of the factor constrained and
#' frees no more than `max.free` items of that factor (Byrne, Shavelson &
#' Muthen, 1989). Whatever is freed stays free in the later steps: an item
#' whose loading is freed also gets a free intercept, and an item whose
#' intercept is freed also gets a free residual variance.
#'
#' Loadings, intercepts and residual variances can be released. Thresholds of
#' ordinal items are not: an item with free thresholds needs its own
#' identification constraints (Wu & Estabrook, 2016), and that is left for
#' later.
#'
#' This is a specification search guided by the data. The released items are
#' a description of where the misfit sits in this sample, not a confirmed
#' finding; the tests of the later steps are conditional on the search.
#'
#' @param x A `psw_invariance` object from [invarianceWiz()].
#' @param auto Search for constraints to release? If `FALSE`, only the items
#'   in `free` are released.
#' @param free Items to release from the start, by name: a list with any of
#'   `load`, `int` and `res`, e.g. `list(int = c("x2", "x3"))`.
#' @param max.free Largest number of items released per factor and parameter
#'   type by the search (default 2).
#' @param min.invariant Smallest number of items per factor that must stay
#'   constrained (default 2).
#' @param criteria,cutoffs As in [invarianceVerdict()]; by default those of `x`.
#'
#' @return An object of class `psw_partial`: a list with data frames `fit`
#'   and `comparison` (as in [invarianceWiz()], plus a column `freed`),
#'   `released` (one row per released constraint, with the score test that
#'   led to it), the list `freed`, a one-paragraph `summary`, and `notes`.
#'
#' @references
#' Byrne, B. M., Shavelson, R. J., & Muthen, B. (1989). Testing for the
#' equivalence of factor covariance and mean structures: The issue of partial
#' measurement invariance. *Psychological Bulletin, 105*(3), 456-466.
#'
#' Steenkamp, J. B. E. M., & Baumgartner, H. (1998). Assessing measurement
#' invariance in cross-national consumer research. *Journal of Consumer
#' Research, 25*(1), 78-90.
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' res <- invarianceWiz(wizdata, f, group = "country")
#' invariancePartial(res)
#' }
#' @export
invariancePartial <- function(x, auto = TRUE, free = list(), max.free = 2L, min.invariant = 2L,
                              criteria = NULL, cutoffs = list()) {
  stopifnot(inherits(x, "psw_invariance"))
  prep <- x$prep
  if (is.null(prep))
    stop("`x` does not carry its prepared data (`x$prep`); run invarianceWiz() again.", call. = FALSE)
  if (is.null(criteria)) criteria <- if (is.null(x$criteria)) c("cfi", "rmsea") else x$criteria
  cut <- utils::modifyList(if (is.null(x$cutoffs)) psw_default_cutoffs else x$cutoffs, as.list(cutoffs))
  it <- prep$items
  ord <- prep$type == "ordinal"
  info <- x$info
  par_arg <- if (ord) info$parameterization else "delta"
  steps <- x$fit$step
  # which list in `freed` each step adds to (thresholds: none, nothing is released there)
  kind_of <- function(s) switch(s, metric = "load", scalar = "int", strict = "res", NA_character_)

  bad <- setdiff(names(free), c("load", "int", "res"))
  if (length(bad)) stop("`free` may only contain `load`, `int` and `res`.", call. = FALSE)
  bad <- setdiff(unlist(free), it$name)
  if (length(bad)) stop("Not an item of the model: ", paste(bad, collapse = ", "), call. = FALSE)
  freed <- list(load = as.character(free$load), int = as.character(free$int), res = as.character(free$res))

  released <- list()
  log <- function(step, item, reason, cu = NULL)
    released[[length(released) + 1L]] <<- data.frame(
      step = psw_step_labels(step, prep$type), parameter = psw_step_target(step), item = item,
      factor = it$factor[match(item, it$name)],
      group = if (is.null(cu)) "" else cu$group, x2 = if (is.null(cu)) NA_real_ else cu$x2,
      df = if (is.null(cu)) NA_integer_ else cu$df, p = if (is.null(cu)) NA_real_ else cu$p,
      reason = reason, stringsAsFactors = FALSE)
  for (s in c("metric", "scalar", "strict"))
    for (i in freed[[kind_of(s)]]) log(s, i, "Chosen by the user")

  fit_step <- function(s)
    psw_inv_task(list(kind = "inv", step = s, model = psw_syntax(prep, s, par_arg, free = freed), prep = prep,
                      estimator = info$estimator, missing = info$missing, parameterization = par_arg))
  judge <- function(comparison)
    invarianceVerdict(structure(list(comparison = comparison), class = "psw_invariance"),
                      criteria = criteria, cutoffs = cut)$comparison
  verdict_of <- function(prev, cur, s)
    judge(psw_ladder_tables(list(prev, cur), c("previous", s), c("previous", s), ord)$comparison)$verdict

  res <- list(configural = fit_step("configural"))
  for (s in steps[-1L]) {
    k <- kind_of(s)
    # what was freed at an earlier step stays free
    carry <- switch(s, scalar = setdiff(freed$load, freed$int), strict = setdiff(freed$int, freed$res), character())
    for (i in carry) log(s, i, if (s == "scalar") "Its loading is free" else "Its intercept is free")
    if (!is.na(k)) freed[[k]] <- c(freed[[k]], carry)

    prev <- res[[length(res)]]
    cur <- fit_step(s)
    if (auto && !is.na(k)) repeat {
      if (!verdict_of(prev, cur, s) %in% c("Not supported", "Ambiguous")) break
      cu <- cur$culprits
      if (is.null(cu)) break
      cu <- cu[!cu$item %in% freed[[k]], , drop = FALSE]
      n_items <- vapply(cu$factor, function(f) sum(it$factor == f), numeric(1))
      n_free <- vapply(cu$factor, function(f) sum(it$factor[match(freed[[k]], it$name)] == f), numeric(1))
      ok <- which(n_items - n_free - 1 >= min.invariant & n_free + 1 <= max.free)
      if (!length(ok)) break
      pick <- cu[ok[1L], , drop = FALSE]
      freed[[k]] <- c(freed[[k]], pick$item)
      log(s, pick$item, "Largest score test at a step that was not supported", pick)
      cur <- fit_step(s)
    }
    res[[s]] <- cur
  }

  freed_at <- vapply(steps, function(s) {
    k <- kind_of(s)
    if (is.na(k)) "" else paste(freed[[k]], collapse = ", ")
  }, "")
  labels <- paste0(psw_step_labels(steps, prep$type), ifelse(nzchar(freed_at), ", partial", ""))
  lad <- psw_ladder_tables(res, steps, labels, ord)
  fit <- cbind(lad$fit, freed = unname(freed_at), stringsAsFactors = FALSE)
  cp <- judge(lad$comparison)
  cp$freed <- unname(freed_at[-1L])

  # how far up the ladder does (partial) invariance go?
  holds <- cp$verdict %in% c("Supported", "Equivalent to previous model")
  reached <- cp$step[seq_len(if (all(holds)) nrow(cp) else which(!holds)[1L] - 1L)]
  listing <- function(k, what) if (length(freed[[k]])) sprintf("%s of %s", what, paste(freed[[k]], collapse = ", "))
  summary <- if ("scalar" %in% reached) {
    fr <- c(listing("load", "loadings"), listing("int", "intercepts"))
    if (length(fr)) paste0(
      "Partial scalar invariance is supported after freeing the ", paste(fr, collapse = " and the "), ". ",
      "With at least two invariant items per factor, latent means are usually considered comparable in a model ",
      "that keeps these parameters free (Byrne, Shavelson & Muthen, 1989). Observed sum or mean scores still ",
      "contain the freed items, so their group differences are not purely differences in the construct.")
    else "Full scalar invariance is supported; no constraint had to be released."
  } else if ("metric" %in% reached) {
    paste0(if (length(freed$load)) sprintf("Partial metric invariance is supported after freeing the %s. ", listing("load", "loadings"))
           else "Metric invariance is supported. ",
           "The intercepts step is not, even with the releases allowed here. Relations among the factors ",
           "(correlations, regressions) can be compared across groups; means should not be.")
  } else {
    sprintf("The %s step is not supported%s, so this model does not justify comparing the groups on these factors.",
            tolower(psw_step_labels(cp$step[which(!holds)[1L]], prep$type)),
            if (cp$step[which(!holds)[1L]] == "thresholds") " (thresholds are not released)"
            else " even with the releases allowed here")
  }

  out <- list(fit = fit, comparison = cp, released = if (length(released)) do.call(rbind, released),
              freed = freed, summary = summary, variant = lad$variant, criteria = criteria, cutoffs = cut,
              means_scalar = if (!is.null(res$scalar)) res$scalar$means else NULL,
              notes = unique(lad$notes))
  class(out) <- "psw_partial"
  out
}

#' Latent mean comparisons across groups
#'
#' Compares factor means across groups from a fitted measurement invariance model,
#' once scalar (or partial scalar) invariance is supported.
#'
#' @param x An object of class `psw_invariance` or `psw_partial`.
#' @param reference Name of the baseline group against which other groups are
#'   compared. Defaults to the first group (`NULL` or `""`).
#' @param alpha Significance level for confidence intervals (default 0.05 for 95% CI).
#'
#' @return An object of class `psw_means`, a data frame containing columns
#'   `factor`, `group`, `is_ref`, `mean`, `se`, `ci_lower`, `ci_upper`, `z`, `p`,
#'   `d_ref`, `d_pool`, with attributes `reference`, `supported`, `model_type`,
#'   `freed_items`, `alpha`, and `notes`.
#'
#' @references
#' Hancock, G. R. (2001). Effect size, power, and sample size determination for
#' structured means modeling and latent growth modeling. *Journal of Educational
#' and Behavioral Statistics*, 26(2), 149-169.
#'
#' Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
#' equivalence of factor covariance and mean structures: The issue of partial
#' measurement invariance. *Psychological Bulletin*, 105(3), 456-466.
#'
#' @export
invarianceMeans <- function(x, reference = NULL, alpha = 0.05) {
  if (!inherits(x, c("psw_invariance", "psw_partial")))
    stop("`x` must be an object of class 'psw_invariance' or 'psw_partial'.", call. = FALSE)

  empty_res <- function(supported = FALSE, note = "", ref = "") {
    df <- data.frame(
      factor = character(), group = character(), is_ref = logical(),
      mean = numeric(), se = numeric(), ci_lower = numeric(), ci_upper = numeric(),
      z = numeric(), p = numeric(), d_ref = numeric(), d_pool = numeric(),
      stringsAsFactors = FALSE
    )
    attr(df, "supported") <- supported
    attr(df, "reference") <- ref
    attr(df, "model_type") <- "none"
    attr(df, "freed_items") <- character()
    attr(df, "alpha") <- alpha
    attr(df, "notes") <- note
    class(df) <- c("psw_means", "data.frame")
    df
  }

  model_type <- "none"
  means_data <- NULL
  freed_items <- character()

  if (inherits(x, "psw_invariance")) {
    cp <- x$comparison
    sc_row <- cp[cp$step == "scalar", , drop = FALSE]
    full_supported <- nrow(sc_row) > 0L && sc_row$verdict[1L] %in% c("Supported", "Equivalent to previous model")

    part <- x$partial
    partial_supported <- FALSE
    if (!is.null(part)) {
      pcp <- part$comparison
      holds <- pcp$verdict %in% c("Supported", "Equivalent to previous model")
      reached <- pcp$step[seq_len(if (all(holds)) nrow(pcp) else which(!holds)[1L] - 1L)]
      partial_supported <- "scalar" %in% reached
    }

    if (full_supported && !is.null(x$means_scalar)) {
      means_data <- x$means_scalar
      model_type <- "full"
    } else if (partial_supported && !is.null(part$means_scalar)) {
      means_data <- part$means_scalar
      model_type <- "partial"
      freed_items <- c(part$freed$load, part$freed$int)
    } else if (!is.null(x$fits$scalar)) {
      prep <- x$prep
      if (is.null(prep)) prep <- list(groups = x$groups$group, n = x$groups$n, fnames = stats::setNames(names(x$info$factors), names(x$info$factors)))
      est <- psw_estimates(x$fits$scalar, prep)
      if (!is.null(est$means)) {
        means_data <- est$means
        model_type <- if (full_supported) "full" else "unsupported"
      }
    }
  } else if (inherits(x, "psw_partial")) {
    pcp <- x$comparison
    holds <- pcp$verdict %in% c("Supported", "Equivalent to previous model")
    reached <- pcp$step[seq_len(if (all(holds)) nrow(pcp) else which(!holds)[1L] - 1L)]
    if ("scalar" %in% reached && !is.null(x$means_scalar)) {
      means_data <- x$means_scalar
      model_type <- "partial"
      freed_items <- c(x$freed$load, x$freed$int)
    }
  }

  if (is.null(means_data) || model_type %in% c("none", "unsupported")) {
    note <- "Latent means cannot be compared: scalar invariance is not supported across groups."
    return(empty_res(supported = FALSE, note = note, ref = if (!is.null(reference)) reference else ""))
  }

  groups <- means_data$groups
  G <- length(groups)
  if (G < 2L) {
    return(empty_res(supported = FALSE, note = "Latent means require at least two groups.", ref = ""))
  }

  if (is.null(reference) || !nzchar(reference)) {
    ref_name <- groups[1L]
  } else {
    if (!reference %in% groups)
      stop(sprintf("Reference group '%s' not found. Available groups: %s", reference, paste(groups, collapse = ", ")), call. = FALSE)
    ref_name <- reference
  }
  ref_idx <- match(ref_name, groups)

  pe <- means_data$pe
  pt <- means_data$pt
  V <- means_data$vcov
  fnames <- means_data$fnames
  N <- means_data$n

  crit <- stats::qnorm(1 - alpha / 2)
  rows <- list()

  for (fid in names(fnames)) {
    fname <- unname(fnames[[fid]])
    m_est <- numeric(G)
    psi_est <- numeric(G)
    psi_est[1L] <- 1.0

    # Group 1: mean is 0 by identification
    m_est[1L] <- 0.0

    # Extract estimates for each group
    for (g in seq_len(G)) {
      r_m <- pe[pe$lhs == fid & pe$op == "~1" & pe$group == g, , drop = FALSE]
      if (nrow(r_m)) m_est[g] <- r_m$est[1L]

      r_psi <- pe[pe$lhs == fid & pe$rhs == fid & pe$op == "~~" & pe$group == g, , drop = FALSE]
      if (nrow(r_psi)) psi_est[g] <- r_psi$est[1L]
    }

    # Variance-covariance matrix of means for this factor
    V_mat <- matrix(0, G, G)
    if (!is.null(V) && !is.null(pt)) {
      for (i in seq(2L, G)) {
        pt_i <- pt[pt$lhs == fid & pt$group == i, , drop = FALSE]
        par_i <- if (nrow(pt_i) && !is.na(pt_i$name[1L])) pt_i$name[1L] else paste0(fid, "~1.g", i)
        if (par_i %in% colnames(V)) {
          for (j in seq(2L, G)) {
            pt_j <- pt[pt$lhs == fid & pt$group == j, , drop = FALSE]
            par_j <- if (nrow(pt_j) && !is.na(pt_j$name[1L])) pt_j$name[1L] else paste0(fid, "~1.g", j)
            if (par_j %in% colnames(V)) {
              V_mat[i, j] <- V[par_i, par_j]
            }
          }
        }
      }
    }

    sd_ref <- sqrt(max(0, psi_est[ref_idx]))
    psi_pool <- sum((N - 1L) * psi_est) / sum(N - 1L)
    sd_pool <- sqrt(max(0, psi_pool))

    for (g in seq_len(G)) {
      is_ref <- (g == ref_idx)
      diff_m <- m_est[g] - m_est[ref_idx]

      if (is_ref) {
        se <- NA_real_
        z <- NA_real_
        p <- NA_real_
        ci_lo <- 0.0
        ci_hi <- 0.0
        d_ref <- 0.0
        d_pool <- 0.0
      } else {
        var_diff <- max(0, V_mat[g, g] + V_mat[ref_idx, ref_idx] - 2 * V_mat[g, ref_idx])
        se <- sqrt(var_diff)
        z <- if (se > 0) diff_m / se else NA_real_
        p <- if (!is.na(z)) 2 * (1 - stats::pnorm(abs(z))) else NA_real_
        ci_lo <- diff_m - crit * se
        ci_hi <- diff_m + crit * se
        d_ref <- if (sd_ref > 0) diff_m / sd_ref else NA_real_
        d_pool <- if (sd_pool > 0) diff_m / sd_pool else NA_real_
      }

      rows[[length(rows) + 1L]] <- data.frame(
        factor = fname,
        group = groups[g],
        is_ref = is_ref,
        mean = diff_m,
        se = se,
        ci_lower = ci_lo,
        ci_upper = ci_hi,
        z = z,
        p = p,
        d_ref = d_ref,
        d_pool = d_pool,
        stringsAsFactors = FALSE
      )
    }
  }

  res <- do.call(rbind, rows)
  rownames(res) <- NULL

  notes <- character()
  if (model_type == "partial") {
    notes <- c(notes, sprintf("Latent means estimated under partial scalar invariance (freed items: %s).",
                              paste(unique(freed_items), collapse = ", ")))
  }

  attr(res, "supported") <- TRUE
  attr(res, "reference") <- ref_name
  attr(res, "model_type") <- model_type
  attr(res, "freed_items") <- freed_items
  attr(res, "alpha") <- alpha
  attr(res, "notes") <- notes
  class(res) <- c("psw_means", "data.frame")
  res
}

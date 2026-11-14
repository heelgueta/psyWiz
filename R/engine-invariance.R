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
#'
#' @return An object of class `psw_invariance`: a list with data frames `fit`
#'   (one row per step), `comparison` (each step against the previous one,
#'   with a verdict), `culprits` (score tests for the constraints of each
#'   step), `loadings`, `groups`, plus `syntax`, `notes` and `info`.
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
                          parallel = "auto", ncores = NULL, keep.fits = FALSE) {
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
    info = list(type = prep$type, estimator = estimator, missing = miss,
                parameterization = parameterization, variant = variant,
                factors = split(prep$items$name, factor(prep$items$factor, levels = unname(prep$fnames))),
                group = group, nthr = stats::setNames(prep$nthr, prep$items$name),
                mode = run$mode, workers = run$workers,
                time = proc.time()[["elapsed"]] - t0,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- lapply(res, function(r) r$fit)
  class(out) <- "psw_invariance"
  invarianceVerdict(out, criteria = criteria, cutoffs = cutoffs)
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

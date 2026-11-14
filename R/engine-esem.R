# engine-esem.R -- CFA against ESEM, and factor correlations under different
# ways of scoring.

psw_esem_task <- function(task) {
  prep <- task$prep
  ids <- prep$items$id
  args <- list(model = task$model, data = prep$data, estimator = task$estimator, missing = task$missing)
  if (task$esem) {
    args$rotation <- task$rotation
    if (task$rotation == "target") args$rotation.args <- list(target = task$target)
  } else {
    args$std.lv <- TRUE
  }
  if (prep$type == "ordinal") args$ordered <- ids
  warns <- character()
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers(do.call(lavaan::cfa, args),
                        warning = function(w) {
                          warns <<- c(warns, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) e)
  r <- list(time = proc.time()[["elapsed"]] - t0, warnings = unique(psw_clean_msg(warns)), error = NULL, fit = NULL)
  if (inherits(fit, "error")) r$error <- psw_clean_msg(conditionMessage(fit)) else r$fit <- fit
  r$converged <- psw_converged(r$fit)
  if (is.null(r$error) && !r$converged) r$error <- "The model did not converge."
  r$fm <- psw_fitmeasures_raw(r$fit)
  if (r$converged) {
    std <- lavaan::lavInspect(r$fit, "std")
    r$lambda <- unclass(std$lambda)
    r$psi <- unclass(std$psi)
    if (isTRUE(task$scores)) {
      r$scores <- tryCatch(suppressWarnings(stats::cor(lavaan::lavPredict(r$fit), use = "pairwise.complete.obs")),
                           error = function(e) NULL)
      # Bartlett scores exist for continuous items only
      r$bartlett <- if (prep$type == "ordinal") NULL else
        tryCatch(suppressWarnings(stats::cor(lavaan::lavPredict(r$fit, method = "Bartlett"),
                                             use = "pairwise.complete.obs")), error = function(e) NULL)
    }
  }
  r
}

# Put ESEM factors in the order (and direction) of the hypothesised factors.
# Each hypothesised factor takes the ESEM factor with the largest mean absolute
# loading on its own items; factors are reflected so those loadings are positive.
psw_esem_align <- function(lambda, psi, target_idx) {
  m <- ncol(lambda)
  score <- vapply(seq_len(m), function(f) colMeans(abs(lambda[target_idx[[f]], , drop = FALSE])), numeric(m))
  perm <- integer(m)
  left <- seq_len(m)
  for (f in order(-apply(score, 2L, max))) {
    pick <- left[which.max(score[left, f])]
    perm[f] <- pick
    left <- setdiff(left, pick)
  }
  lambda <- lambda[, perm, drop = FALSE]
  psi <- psi[perm, perm, drop = FALSE]
  sgn <- vapply(seq_len(m), function(f) if (sum(lambda[target_idx[[f]], f]) < 0) -1 else 1, numeric(1))
  lambda <- sweep(lambda, 2L, sgn, `*`)
  psi <- psi * outer(sgn, sgn)
  list(lambda = lambda, psi = psi)
}

#' CFA against ESEM, and factor correlations by scoring method
#'
#' Fits the hypothesised correlated-factors model twice: as a CFA (every item
#' loads on its own factor only) and as an exploratory structural equation
#' model (ESEM), where every item may load on every factor and the solution
#' is rotated towards the hypothesised pattern. Comparing the two shows
#' whether forcing cross-loadings to zero distorts the picture, which
#' typically inflates factor correlations.
#'
#' The second part answers a different everyday question: how much does the
#' correlation between two constructs depend on how they are scored? It
#' lists every factor correlation as obtained from composite (mean) scores,
#' composites corrected for unreliability (with alpha), the CFA, CFA factor
#' scores (regression method, and Bartlett's method for continuous items),
#' and the ESEM.
#'
#' @inheritParams invarianceWiz
#' @param factors A named list with at least two factors.
#' @param rotation `"target"` (default: rotate towards the hypothesised
#'   pattern, cross-loadings towards zero) or `"geomin"` (no hypothesis used
#'   in the rotation).
#' @return An object of class `psw_esem`: a list with data frames `fit`
#'   (CFA and ESEM), `comparison` (chi-square difference test), `loadings`
#'   (one row per item and factor, CFA and ESEM standardised loadings),
#'   `correlations` (one row per factor pair, one column per method), plus
#'   `notes` and `info`.
#' @references
#' Asparouhov, T., & Muthen, B. (2009). Exploratory structural equation
#' modeling. *Structural Equation Modeling, 16*(3), 397-438.
#'
#' Marsh, H. W., Morin, A. J. S., Parker, P. D., & Kaur, G. (2014).
#' Exploratory structural equation modeling: An integration of the best
#' features of exploratory and confirmatory factor analysis. *Annual Review
#' of Clinical Psychology, 10*, 85-110.
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' esemWiz(wizdata, f)
#' }
#' @export
esemWiz <- function(data, factors, type = c("auto", "ordinal", "continuous"), estimator = "auto",
                    missing = c("auto", "listwise", "fiml", "pairwise"),
                    rotation = c("target", "geomin"),
                    parallel = "auto", ncores = NULL, keep.fits = FALSE) {
  type <- match.arg(type)
  missing <- match.arg(missing)
  rotation <- match.arg(rotation)
  prep <- psw_prepare(data, factors, type = type)
  m <- length(prep$fnames)
  if (m < 2L) stop("ESEM needs at least two factors; with one factor it is the same model as the CFA.", call. = FALSE)
  ord <- prep$type == "ordinal"
  notes <- prep$notes
  it <- prep$items

  estimator <- toupper(estimator)
  if (estimator == "AUTO") estimator <- if (ord) "WLSMV" else "MLR"
  allowed <- if (ord) c("WLSMV", "ULSMV", "DWLS", "ULS") else psw_ml_estimators
  if (!estimator %in% allowed)
    stop(sprintf("Estimator '%s' is not available for %s items. Choose one of: %s.",
                 estimator, prep$type, paste(allowed, collapse = ", ")), call. = FALSE)
  miss <- psw_resolve_missing(missing, estimator, ord, prep$has_missing)
  if (ord && miss == "fiml") stop("FIML is not available for ordinal items.", call. = FALSE)
  if (!ord && miss == "pairwise") stop("Pairwise deletion is only available for ordinal items.", call. = FALSE)

  fids <- names(prep$fnames)
  target_idx <- lapply(fids, function(f) which(it$fid == f))
  target <- matrix(0, nrow(it), m)
  for (f in seq_len(m)) target[target_idx[[f]], f] <- NA_real_
  cfa_syntax <- psw_syntax(prep)
  esem_syntax <- paste0(paste0('efa("b")*', fids, collapse = " + "), " =~ ", paste(it$id, collapse = " + "))
  common <- list(kind = "esem", prep = prep, estimator = estimator, missing = miss)
  tasks <- list(c(common, list(model = cfa_syntax, esem = FALSE, scores = TRUE)),
                c(common, list(model = esem_syntax, esem = TRUE, rotation = rotation, target = target)))

  t0 <- proc.time()[["elapsed"]]
  run <- psw_run_tasks(tasks, parallel = parallel, ncores = ncores)
  res <- stats::setNames(run$results, c("cfa", "esem"))
  labels <- c("CFA", sprintf("ESEM (%s rotation)", rotation))
  # the ladder helper compares each model with the previous one, and the
  # chi-square difference test needs the less restricted model first
  lad <- psw_ladder_tables(res[c("esem", "cfa")], c("esem", "cfa"), rev(labels), ord)
  notes <- c(notes, lad$notes)
  fit <- lad$fit[2:1, , drop = FALSE]
  rownames(fit) <- NULL

  # ---- loadings: every item on every factor ----
  fn <- unname(prep$fnames)
  grid <- expand.grid(i = seq_len(nrow(it)), f = seq_len(m))
  loadings <- data.frame(item = it$name[grid$i], own = it$factor[grid$i], factor = fn[grid$f],
                         target = it$fid[grid$i] == fids[grid$f], cfa = NA_real_, esem = NA_real_,
                         stringsAsFactors = FALSE)
  esem_psi <- NULL
  if (res$cfa$converged) {
    L <- res$cfa$lambda
    loadings$cfa <- ifelse(loadings$target, L[cbind(grid$i, grid$f)], NA_real_)
  }
  if (res$esem$converged) {
    al <- psw_esem_align(res$esem$lambda, res$esem$psi, target_idx)
    loadings$esem <- al$lambda[cbind(grid$i, grid$f)]
    esem_psi <- al$psi
  }

  # ---- factor correlations by method ----
  num <- if (ord) psw_prepare(data, factors, type = "continuous") else prep
  comp <- vapply(fids, function(f) rowMeans(num$data[num$items$id[num$items$fid == f]], na.rm = TRUE),
                 numeric(nrow(num$data)))
  rc <- stats::cor(comp, use = "pairwise.complete.obs")
  alpha <- vapply(fids, function(f) {
    X <- num$data[num$items$id[num$items$fid == f]]
    psw_alpha(stats::cov(X, use = "complete.obs"))
  }, numeric(1))
  pairs <- which(upper.tri(diag(m)), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  pick <- function(M) if (is.null(M)) rep(NA_real_, nrow(pairs)) else M[pairs]
  correlations <- data.frame(
    factor1 = fn[pairs[, 1L]], factor2 = fn[pairs[, 2L]],
    composite = pick(rc),
    corrected = pick(rc / sqrt(outer(alpha, alpha))),
    cfa = pick(res$cfa$psi), scores = pick(res$cfa$scores), bartlett = pick(res$cfa$bartlett),
    esem = pick(esem_psi),
    stringsAsFactors = FALSE)
  if (!is.null(lad$comparison) && isTRUE(lad$variant[["chisq"]] == "scaled"))
    notes <- c(notes, "The chi-square difference test is a scaled test; it is not the difference between the two chi-square values shown.")
  if (any(alpha <= 0, na.rm = TRUE))
    notes <- c(notes, "A scale has alpha at or below zero, so the corrected composite correlation is not meaningful.")

  out <- list(
    fit = fit, comparison = lad$comparison, loadings = loadings, correlations = correlations,
    reliability = data.frame(factor = fn, alpha = unname(alpha), stringsAsFactors = FALSE),
    syntax = list(cfa = psw_syntax_display(cfa_syntax, prep), esem = psw_syntax_display(esem_syntax, prep)),
    notes = unique(notes),
    info = list(type = prep$type, estimator = estimator, missing = miss, rotation = rotation,
                variant = lad$variant, factors = split(it$name, factor(it$factor, levels = fn)),
                n = nrow(prep$data), mode = run$mode, workers = run$workers,
                time = proc.time()[["elapsed"]] - t0,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- lapply(res, function(r) r$fit)
  class(out) <- "psw_esem"
  out
}

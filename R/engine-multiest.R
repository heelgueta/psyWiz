# engine-multiest.R -- the same CFA under several estimators, side by side.

# One estimator; runs in the main process or in a worker.
psw_est_task <- function(task) {
  r <- psw_fit(task$model, task$prep, task$estimator, task$missing)
  r$estimator <- task$estimator
  r$converged <- psw_converged(r$fit)
  if (is.null(r$error) && !r$converged) r$error <- "The model did not converge."
  r$fm <- psw_fitmeasures_raw(r$fit)
  est <- psw_estimates(r$fit, task$prep)
  r$loadings <- est$loadings
  r$cors <- est$cors
  r$reliability <- psw_reliability(est$loadings)
  r
}

#' Confirmatory factor analysis under several estimators
#'
#' Fits one single-group CFA (one factor or several correlated factors) with
#' each of the requested estimators and returns the fit indices, loadings,
#' factor correlations and reliabilities side by side.
#'
#' Maximum-likelihood estimators (`"ML"`, `"MLR"`, `"MLM"`) always treat the
#' items as continuous. Least-squares estimators (`"ULS"`, `"ULSMV"`,
#' `"DWLS"`, `"WLSMV"`, `"WLS"`) treat them as ordinal (polychoric
#' correlations) when the items are ordinal, and as continuous otherwise. This
#' is the comparison applied researchers usually want: "what changes if I
#' treat my Likert items as continuous?"
#'
#' @inheritParams invarianceWiz
#' @param estimators Character vector of estimators, or `"MULTI"` for all of
#'   `"ML"`, `"MLR"`, `"MLM"`, `"ULS"`, `"ULSMV"`, `"DWLS"`, `"WLSMV"`, `"WLS"`.
#' @param missing `"auto"`, `"listwise"`, `"fiml"` or `"pairwise"`. A method an
#'   estimator cannot use is replaced by listwise deletion for that estimator
#'   (with a note). `"auto"` uses FIML for ML/MLR and pairwise deletion for
#'   ordinal treatment when there are missing values.
#' @param parallel `"auto"`, `TRUE` or `FALSE`: fit the estimators on several
#'   cores? `"auto"` does so only when the first model is slow.
#'
#' @return An object of class `psw_multi`: a list with data frames `fit` (one
#'   row per estimator), `loadings`, `cors` and `reliability` (long format,
#'   one block per estimator), plus `syntax`, `notes` and `info`.
#'
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' cfaWiz(wizdata, f, estimators = c("MLR", "WLSMV"))
#' }
#' @export
cfaWiz <- function(data, factors, estimators = c("MLR", "ULSMV", "WLSMV"),
                   type = c("auto", "ordinal", "continuous"),
                   missing = c("auto", "listwise", "fiml", "pairwise"),
                   parallel = "auto", ncores = NULL, keep.fits = FALSE) {
  type <- match.arg(type)
  missing <- match.arg(missing)
  estimators <- toupper(as.character(estimators))
  if (any(estimators %in% c("MULTI", "ALL"))) estimators <- psw_all_estimators
  estimators <- unique(estimators)
  bad <- setdiff(estimators, psw_all_estimators)
  if (length(bad))
    stop("Unknown estimator(s): ", paste(bad, collapse = ", "), ". Available: ",
         paste(psw_all_estimators, collapse = ", "), ".", call. = FALSE)
  if (!length(estimators)) stop("Choose at least one estimator.", call. = FALSE)

  base <- psw_prepare(data, factors, type = type)
  preps <- list(continuous = if (base$type == "continuous") base
                             else psw_prepare(data, factors, type = "continuous"))
  if (base$type == "ordinal") preps$ordinal <- base
  notes <- base$notes
  syntax <- psw_syntax(base)

  treated <- ifelse(base$type == "ordinal" & estimators %in% psw_ls_estimators, "ordinal", "continuous")
  miss <- character(length(estimators))
  for (k in seq_along(estimators)) {
    e <- estimators[k]
    ordk <- treated[k] == "ordinal"
    m <- psw_resolve_missing(missing, e, ordk, preps[[treated[k]]]$has_missing)
    valid <- (m == "listwise") || (m == "fiml" && e %in% c("ML", "MLR")) || (m == "pairwise" && ordk)
    if (!valid) {
      notes <- c(notes, sprintf("%s cannot use '%s' for missing values with items treated as %s; listwise deletion was used.",
                                e, m, treated[k]))
      m <- "listwise"
    }
    miss[k] <- m
  }

  tasks <- lapply(seq_along(estimators), function(k)
    list(kind = "est", estimator = estimators[k], model = syntax, prep = preps[[treated[k]]], missing = miss[k]))

  t0 <- proc.time()[["elapsed"]]
  run <- psw_run_tasks(tasks, parallel = parallel, ncores = ncores)
  res <- run$results

  fit <- do.call(rbind, lapply(seq_along(estimators), function(k) {
    r <- res[[k]]
    v <- psw_choose_variant(list(r$fm), treated[k] == "ordinal")
    data.frame(estimator = estimators[k], treated = treated[k], missing = miss[k],
               converged = r$converged, t(psw_measures(r$fm, v)),
               v.chisq = v[["chisq"]], v.cfi = v[["cfi"]], v.rmsea = v[["rmsea"]],
               time = r$time, stringsAsFactors = FALSE)
  }))
  rownames(fit) <- NULL
  if (any(estimators %in% c("ULS", "DWLS") & fit$converged))
    notes <- c(notes, paste0("ULS and DWLS without a robust correction give no valid p value for the ",
                             "chi-square test; ULSMV and WLSMV are the corrected versions."))
  if (length(unique(fit$n[fit$converged])) > 1L)
    notes <- c(notes, "The estimators did not all use the same number of cases (see N), because they handle missing values differently.")
  for (k in seq_along(estimators)) {
    r <- res[[k]]
    if (!is.null(r$error)) notes <- c(notes, sprintf("%s: %s", estimators[k], r$error))
    if (length(r$warnings)) notes <- c(notes, sprintf("%s: %s", estimators[k], r$warnings))
  }
  stack <- function(what) {
    out <- do.call(rbind, lapply(seq_along(estimators), function(k) {
      d <- res[[k]][[what]]
      if (is.null(d) || nrow(d) == 0L) return(NULL)
      d$group <- NULL
      cbind(estimator = estimators[k], d, stringsAsFactors = FALSE)
    }))
    if (!is.null(out)) rownames(out) <- NULL
    out
  }

  out <- list(
    fit = fit, loadings = stack("loadings"), cors = stack("cors"), reliability = stack("reliability"),
    syntax = psw_syntax_display(syntax, base), notes = unique(notes),
    info = list(type = base$type, estimators = estimators,
                factors = split(base$items$name, factor(base$items$factor, levels = unname(base$fnames))),
                n = nrow(base$data), mode = run$mode, workers = run$workers,
                time = proc.time()[["elapsed"]] - t0,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- stats::setNames(lapply(res, function(r) r$fit), estimators)
  class(out) <- "psw_multi"
  out
}

# Long estimate table -> wide matrix (rows = `rows`, columns = estimators).
psw_wide <- function(long, row_cols, value, estimators) {
  if (is.null(long) || nrow(long) == 0L) return(NULL)
  key <- do.call(paste, c(long[row_cols], sep = "\r"))
  rows <- unique(key)
  m <- matrix(NA_real_, length(rows), length(estimators), dimnames = list(NULL, estimators))
  m[cbind(match(key, rows), match(long$estimator, estimators))] <- long[[value]]
  first <- match(rows, key)
  cbind(long[first, row_cols, drop = FALSE], as.data.frame(m, check.names = FALSE), row.names = NULL)
}

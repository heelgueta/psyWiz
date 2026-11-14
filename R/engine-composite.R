# engine-composite.R -- is ONE composite an adequate summary of a set of
# formative indicators?
#
# Reflective indicators can be checked against each other (they should
# correlate because they share a cause). Formative indicators need not
# correlate, so internal-consistency logic says nothing about them, and the
# usual advice is to look at how well the composite predicts something. That
# cannot fail: any weighted sum predicts something.
#
# What CAN fail is the claim that the indicators act on other variables
# *through a single composite*. If that is true, each indicator's effects on
# different outcomes must be proportional (Franke, Preacher & Rigdon, 2008):
# the p x q matrix of regression coefficients of q outcomes on p indicators
# has rank one. That is exactly the reduced-rank regression model of
# econometrics (Anderson, 1951), and the rank can be tested with the
# canonical correlations between indicators and outcomes (Bartlett, 1938).
#
# compositeWiz() therefore compares three ways of predicting the outcomes:
#   free      every outcome has its own weights for the indicators
#   composite one set of weights, estimated (rank-one reduced-rank regression)
#   unit      one set of weights, all equal (the plain sum most people compute)
# in sample and, because the free model always fits best in sample, in
# cross-validation.

psw_r2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

# Rank-one weights: the first canonical variate of X with respect to Y.
psw_rank1_weights <- function(X, Y) {
  cc <- stats::cancor(X, Y, xcenter = TRUE, ycenter = TRUE)
  list(w = cc$xcoef[, 1L], rho = cc$cor)
}

psw_composite_predict <- function(Xtr, Ytr, Xte, method) {
  mx <- colMeans(Xtr)
  sx <- apply(Xtr, 2L, stats::sd)
  Ztr <- scale(Xtr, mx, sx)
  Zte <- scale(Xte, mx, sx)
  if (method == "free") {
    B <- stats::coef(stats::lm.fit(cbind(1, Ztr), Ytr))
    return(cbind(1, Zte) %*% B)
  }
  w <- if (method == "unit") rep(1, ncol(Ztr)) else psw_rank1_weights(Ztr, Ytr)$w
  ctr <- drop(Ztr %*% w)
  cte <- drop(Zte %*% w)
  B <- stats::coef(stats::lm.fit(cbind(1, ctr), Ytr))
  cbind(1, cte) %*% B
}

#' Is one composite enough? A check for formative indicators
#'
#' Tests whether a set of formative indicators relates to two or more
#' outcome variables through a single weighted composite, and how much is
#' lost by using equal weights.
#'
#' Three predictions of the outcomes are compared: from the indicators with
#' separate weights for every outcome ("free"), from one composite with
#' estimated weights (rank-one reduced-rank regression), and from the
#' unit-weighted sum of the standardised indicators. If the indicators form
#' one construct, the single composite should predict every outcome about as
#' well as the free model, and the rank test should not reject.
#'
#' Because the free model always wins in sample, the comparison is repeated
#' with k-fold cross-validation; negative out-of-sample R-squared means
#' "worse than predicting the mean".
#'
#' Indicators should be scored in the same conceptual direction before
#' computing the unit-weighted sum.
#'
#' @param data A data frame.
#' @param indicators Names of at least two numeric formative indicators.
#' @param outcomes Names of at least two numeric outcome variables. With a
#'   single outcome the single-composite hypothesis cannot fail.
#' @param folds Number of cross-validation folds.
#' @param repeats Number of times the cross-validation is repeated.
#' @param seed Random seed for the fold assignment.
#' @return An object of class `psw_composite`: a list with `test` (rank
#'   test), `r2` (per outcome: in-sample and cross-validated R-squared under
#'   the three models), `weights` (per indicator: composite weight, its
#'   correlation with the composite, VIF), `summary`, `notes` and `info`.
#' @references
#' Anderson, T. W. (1951). Estimating linear restrictions on regression
#' coefficients for multivariate normal distributions. *Annals of
#' Mathematical Statistics, 22*(3), 327-351.
#'
#' Franke, G. R., Preacher, K. J., & Rigdon, E. E. (2008). Proportional
#' structural effects of formative indicators. *Journal of Business
#' Research, 61*(12), 1229-1237.
#' @export
compositeWiz <- function(data, indicators, outcomes, folds = 10L, repeats = 5L, seed = 1L) {
  indicators <- unique(as.character(indicators))
  outcomes <- unique(as.character(outcomes))
  miss <- setdiff(c(indicators, outcomes), names(data))
  if (length(miss)) stop("Variable(s) not found in the data: ", paste(miss, collapse = ", "), call. = FALSE)
  if (length(indicators) < 2L) stop("At least two indicators are needed.", call. = FALSE)
  if (length(outcomes) < 2L)
    stop("At least two outcomes are needed: with one outcome a single composite fits by construction.", call. = FALSE)
  if (length(intersect(indicators, outcomes))) stop("A variable cannot be both an indicator and an outcome.", call. = FALSE)
  d <- as.data.frame(lapply(data[c(indicators, outcomes)], psw_as_numeric))
  names(d) <- c(indicators, outcomes)
  cc <- stats::complete.cases(d)
  notes <- character()
  if (any(!cc)) notes <- c(notes, sprintf("%d case(s) with missing values were dropped.", sum(!cc)))
  d <- d[cc, , drop = FALSE]
  n <- nrow(d)
  p <- length(indicators)
  q <- length(outcomes)
  if (n < 5L * (p + q)) notes <- c(notes, "The sample is small for the number of variables; treat the results as rough.")
  if (n <= p + q + 2L) stop("Too few complete cases.", call. = FALSE)
  X <- as.matrix(d[indicators])
  Y <- as.matrix(d[outcomes])
  if (any(apply(cbind(X, Y), 2L, stats::sd) == 0)) stop("A variable has no variance.", call. = FALSE)
  Z <- scale(X)

  # ---- rank test (Bartlett): are canonical correlations 2..min(p, q) zero? ----
  r1 <- psw_rank1_weights(Z, Y)
  rho <- r1$rho
  rest <- rho[-1L]
  stat <- -(n - 1 - (p + q + 1) / 2) * sum(log(1 - rest^2))
  df <- (p - 1) * (q - 1)
  test <- data.frame(chisq = stat, df = df, p = stats::pchisq(stat, df, lower.tail = FALSE),
                     rho1 = rho[1L], rho2 = if (length(rho) > 1L) rho[2L] else NA_real_)

  # ---- weights ----
  w <- r1$w
  comp <- drop(Z %*% w)
  if (sum(stats::cor(Z, comp)) < 0) {
    w <- -w
    comp <- -comp
  }
  w <- w / stats::sd(comp)
  comp <- comp / stats::sd(comp)
  vif <- diag(solve(stats::cor(X)))
  weights <- data.frame(indicator = indicators, weight = unname(w),
                        loading = as.numeric(stats::cor(Z, comp)), vif = unname(vif),
                        stringsAsFactors = FALSE)
  if (any(weights$weight * weights$loading < 0))
    notes <- c(notes, paste0("An indicator's weight and its correlation with the composite have opposite signs ",
                             "(suppression, usually from collinearity); its weight should not be interpreted."))

  # ---- R-squared in sample and cross-validated ----
  methods <- c("free", "composite", "unit")
  ins <- vapply(methods, function(m) {
    pred <- psw_composite_predict(X, Y, X, m)
    vapply(seq_len(q), function(k) psw_r2(Y[, k], pred[, k]), numeric(1))
  }, numeric(q))
  folds <- max(2L, min(as.integer(folds), n %/% (p + 3L)))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  cv <- array(NA_real_, c(q, length(methods), repeats))
  for (r in seq_len(repeats)) {
    fold <- sample(rep(seq_len(folds), length.out = n))
    pred <- array(NA_real_, c(n, q, length(methods)))
    for (f in seq_len(folds)) {
      te <- fold == f
      for (m in seq_along(methods))
        pred[te, , m] <- psw_composite_predict(X[!te, , drop = FALSE], Y[!te, , drop = FALSE],
                                               X[te, , drop = FALSE], methods[m])
    }
    for (m in seq_along(methods))
      cv[, m, r] <- vapply(seq_len(q), function(k) psw_r2(Y[, k], pred[, k, m]), numeric(1))
  }
  cvm <- apply(cv, c(1L, 2L), mean)
  r2 <- data.frame(outcome = outcomes,
                   free = ins[, "free"], composite = ins[, "composite"], unit = ins[, "unit"],
                   cv.free = cvm[, 1L], cv.composite = cvm[, 2L], cv.unit = cvm[, 3L],
                   stringsAsFactors = FALSE, row.names = NULL)
  r2$retained <- r2$composite / r2$free

  best <- methods[apply(cvm, 1L, which.max)]
  summary <- c(
    if (test$p < 0.05)
      sprintf("The single-composite hypothesis is rejected (chi-square(%d) = %.2f, p %s): the indicators relate to the outcomes in more than one way.",
              df, stat, if (test$p < 0.001) "< .001" else sprintf("= %.3f", test$p))
    else
      sprintf("The single-composite hypothesis is not rejected (chi-square(%d) = %.2f, p = %.3f): the indicators' effects are compatible with one composite.",
              df, stat, test$p),
    sprintf("One estimated composite retains %.0f%% to %.0f%% of the variance the indicators explain in each outcome.",
            100 * min(r2$retained), 100 * max(r2$retained)),
    sprintf("Out of sample, the best predictor is: %s.",
            paste(sprintf("%s for %s", c(free = "separate weights", composite = "the estimated composite",
                                         unit = "the unit-weighted sum")[best], outcomes), collapse = "; ")))
  out <- list(test = test, r2 = r2, weights = weights, summary = summary, notes = notes,
              info = list(n = n, indicators = indicators, outcomes = outcomes, folds = folds, repeats = repeats))
  class(out) <- "psw_composite"
  out
}

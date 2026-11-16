# engine-bifactor.R -- bifactor CFA: every item loads on one general factor and
# on one specific factor, and all factors are uncorrelated.
#
# The question a bifactor model answers in scale work: is the questionnaire
# essentially one thing (report a total score), or do the subscales carry
# reliable information of their own (report subscale scores)? The loadings
# alone do not answer it; the indices below do (Rodriguez, Reise & Haviland,
# 2016). With standardized loadings g_i (general) and s_i (specific) and
# unique variances u_i = 1 - g_i^2 - s_i^2:
#
#   omega (total score)   [(sum g)^2 + sum_k (sum_{i in k} s)^2] / V
#   omega hierarchical    (sum g)^2 / V             V = numerator of omega + sum u
#   omega (subscale k)    [(sum_k g)^2 + (sum_k s)^2] / V_k
#   omega hier. subscale  (sum_k s)^2 / V_k         share of the subscale score that is NOT general
#   ECV                   sum g^2 / (sum g^2 + sum s^2)      explained common variance
#   item ECV              g_i^2 / (g_i^2 + s_i^2)
#   PUC                   share of item pairs from different specific factors, whose
#                         correlation is due to the general factor only
#
# With ordinal items these describe the latent responses, as everywhere else.

#' Bifactor CFA
#'
#' Fits a bifactor model: one general factor on which every item loads, plus
#' one specific factor per group of items, all uncorrelated. Reports the fit
#' next to that of the correlated-factors model, the loadings on both
#' factors, and the indices that say how much of the scores is general:
#' omega, omega hierarchical, explained common variance.
#'
#' @inheritParams invarianceWiz
#' @param factors A named list with at least two specific factors of three or
#'   more items each.
#' @param general Name of the general factor.
#' @return An object of class `psw_bifactor`: a list with data frames `fit`
#'   (bifactor and correlated-factors models), `loadings` (one row per item:
#'   `general`, `own`, `unique`, `iecv`), `indices`, plus `syntax`, `notes`
#'   and `info`.
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. *Psychological
#' Methods, 21*(2), 137-150.
#' @examples
#' \donttest{
#' data(wizdata)
#' # needs three or more items per specific factor; see plot(x, type = "antenna")
#' }
#' @export
bifactorWiz <- function(data, factors, general = "General", type = c("auto", "ordinal", "continuous"),
                        estimator = "auto", missing = c("auto", "listwise", "fiml", "pairwise"),
                        keep.fits = FALSE) {
  type <- match.arg(type)
  missing <- match.arg(missing)
  prep <- psw_prepare(data, factors, type = type)
  it <- prep$items
  fids <- names(prep$fnames)
  fn <- unname(prep$fnames)
  if (length(fids) < 2L) stop("A bifactor model needs at least two specific factors.", call. = FALSE)
  if (any(table(it$fid) < 3L))
    stop("Every specific factor needs at least three items (with two, its loadings are not identified).", call. = FALSE)
  ord <- prep$type == "ordinal"
  est <- psw_resolve_estimation(estimator, missing, prep)
  notes <- prep$notes

  pairs <- utils::combn(fids, 2L)
  syntax <- paste(c(
    paste0("g =~ ", paste(it$id, collapse = " + ")),
    vapply(fids, function(f) paste0(f, " =~ ", paste(it$id[it$fid == f], collapse = " + ")), ""),
    paste0("g ~~ 0*", fids), sprintf("%s ~~ 0*%s", pairs[1L, ], pairs[2L, ])), collapse = "\n")
  r <- psw_fit(syntax, prep, est$estimator, est$missing)
  if (!is.null(r$error)) stop(r$error, call. = FALSE)
  if (!psw_converged(r$fit))
    stop("The bifactor model did not converge. This is common when a specific factor has little variance of its own.", call. = FALSE)
  notes <- c(notes, r$warnings)
  # the usual alternative, for comparison
  cf <- psw_fit(psw_syntax(prep), prep, est$estimator, est$missing)
  fms <- list(psw_fitmeasures_raw(r$fit), psw_fitmeasures_raw(cf$fit))
  variant <- psw_choose_variant(fms, ord)
  fit <- data.frame(label = c("Bifactor", "Correlated factors"),
                    do.call(rbind, lapply(fms, function(fm) t(psw_measures(fm, variant)))), stringsAsFactors = FALSE)
  if (is.null(fms[[2L]])) notes <- c(notes, "The correlated-factors model did not converge.")

  pe <- lavaan::parameterEstimates(r$fit, standardized = TRUE, ci = FALSE)
  pe <- pe[pe$op == "=~", , drop = FALSE]
  g <- pe$std.all[pe$lhs == "g"][match(it$id, pe$rhs[pe$lhs == "g"])]
  s <- pe$std.all[pe$lhs != "g"][match(it$id, pe$rhs[pe$lhs != "g"])]
  loadings <- data.frame(specific = it$factor, item = it$name, general = g, own = s,
                         unique = 1 - g^2 - s^2, iecv = g^2 / (g^2 + s^2), stringsAsFactors = FALSE)
  if (any(loadings$unique < 0, na.rm = TRUE))
    notes <- c(notes, "An item has a negative unique variance (a Heywood case): the solution is not trustworthy.")

  by <- function(v) vapply(fn, function(f) sum(v[it$factor == f]), numeric(1))
  sg <- by(g); ss <- by(s); su <- by(loadings$unique)
  V <- sum(g)^2 + sum(ss^2) + sum(su)
  Vk <- sg^2 + ss^2 + su
  indices <- rbind(
    data.frame(scale = "Total score", k = nrow(it), omega = (sum(g)^2 + sum(ss^2)) / V, omegaH = sum(g)^2 / V,
               ecv = sum(g^2) / sum(g^2 + s^2), stringsAsFactors = FALSE),
    data.frame(scale = fn, k = as.integer(table(factor(it$factor, levels = fn))), omega = unname((sg^2 + ss^2) / Vk),
               omegaH = unname(ss^2 / Vk), ecv = unname(by(s^2) / (by(g^2) + by(s^2))), stringsAsFactors = FALSE))
  k <- as.numeric(table(it$fid))
  puc <- 1 - sum(k * (k - 1) / 2) / (nrow(it) * (nrow(it) - 1) / 2)

  out <- list(
    fit = fit, loadings = loadings, indices = indices,
    syntax = gsub("(^|\n)g ", paste0("\\1", general, " "), psw_syntax_display(syntax, prep)), notes = unique(notes),
    info = list(type = prep$type, estimator = est$estimator, missing = est$missing, variant = variant,
                general = general, factors = split(it$name, factor(it$factor, levels = fn)), puc = puc,
                n = nrow(prep$data), mode = "sequential", workers = 0L, time = r$time + cf$time,
                lavaan = as.character(utils::packageVersion("lavaan"))))
  if (keep.fits) out$fits <- list(bifactor = r$fit, correlated = cf$fit)
  class(out) <- "psw_bifactor"
  out
}

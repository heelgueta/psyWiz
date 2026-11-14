# engine-fit.R -- fitting one model and pulling plain tables out of it.

psw_ml_estimators <- c("ML", "MLR", "MLM")
psw_ls_estimators <- c("ULS", "ULSMV", "DWLS", "WLSMV", "WLS")
psw_all_estimators <- c(psw_ml_estimators, psw_ls_estimators)

# Resolve missing = "auto" for one estimator / item treatment.
psw_resolve_missing <- function(missing, estimator, ordinal, has_missing) {
  if (missing != "auto") return(missing)
  if (!has_missing) return("listwise")
  if (ordinal) return("pairwise")
  if (estimator %in% c("ML", "MLR")) return("fiml")
  "listwise"
}

#' Fit one lavaan model, capturing warnings and errors (internal)
#' @return list(fit, error, warnings, time)
#' @keywords internal
#' @noRd
psw_fit <- function(model, prep, estimator, missing = "listwise", parameterization = "delta") {
  ids <- prep$items$id
  args <- list(model = model, data = prep$data, std.lv = TRUE,
               estimator = estimator, missing = missing)
  if (length(prep$groups) >= 2L) {
    args$group <- psw_group_col
    args$group.label <- prep$groups
    args$meanstructure <- TRUE
  }
  if (prep$type == "ordinal") {
    args$ordered <- ids
    args$parameterization <- parameterization
  }

  warns <- character()
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers(do.call(lavaan::cfa, args),
                        warning = function(w) {
                          warns <<- c(warns, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) e)
  time <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error"))
    return(list(fit = NULL, error = psw_clean_msg(conditionMessage(fit)), warnings = warns, time = time))
  list(fit = fit, error = NULL, warnings = unique(psw_clean_msg(warns)), time = time)
}

# lavaan messages come with a "lavaan->fn():" prefix and hard line breaks.
psw_clean_msg <- function(x) {
  x <- gsub("\\s+", " ", x)
  x <- sub("^\\s*lavaan(->[A-Za-z0-9_.]+\\(\\))?:\\s*", "", x)
  x <- sub("^lavaan (WARNING|ERROR|NOTE):\\s*", "", x)
  trimws(x)
}

psw_converged <- function(fit) {
  !is.null(fit) && isTRUE(tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE))
}

# All fit measures lavaan can compute, as a plain named vector (NULL if not converged).
psw_fitmeasures_raw <- function(fit) {
  if (!psw_converged(fit)) return(NULL)
  fm <- tryCatch(suppressWarnings(lavaan::fitMeasures(fit)), error = function(e) NULL)
  if (is.null(fm)) return(NULL)
  fm <- unlist(unclass(fm))
  fm[["ntotal"]] <- tryCatch(lavaan::lavInspect(fit, "ntotal"), error = function(e) NA_real_)
  fm
}

psw_variant_suffixes <- c(robust = ".robust", scaled = ".scaled", standard = "")

#' Choose which version of each fit index to report (internal)
#'
#' lavaan offers up to three versions of the chi-square based indices:
#' standard, scaled and robust. psyWiz reports one version per index and uses
#' the same version for every model in a table, so that differences between
#' models never mix versions. Preference: robust > scaled > standard for
#' continuous items; scaled > standard for ordinal items (lavaan's robust CFI
#' and RMSEA for ordinal data are experimental and often unavailable).
#'
#' @param fms List of raw fit-measure vectors (NULL entries are ignored).
#' @param ordinal Were the items treated as ordinal?
#' @return Named character vector (chisq, cfi, rmsea) of variant names.
#' @keywords internal
#' @noRd
psw_choose_variant <- function(fms, ordinal) {
  fms <- Filter(Negate(is.null), fms)
  out <- c(chisq = NA_character_, cfi = NA_character_, rmsea = NA_character_)
  if (!length(fms)) return(out)
  ok <- function(base, variant) {
    nm <- paste0(base, psw_variant_suffixes[[variant]])
    all(vapply(fms, function(fm) nm %in% names(fm) && is.finite(fm[[nm]]), logical(1)))
  }
  first <- function(base, variants) {
    for (v in variants) if (ok(base, v)) return(v)
    NA_character_
  }
  idx <- if (ordinal) c("scaled", "standard") else c("robust", "scaled", "standard")
  out[["chisq"]] <- first("chisq", c("scaled", "standard"))
  out[["cfi"]] <- first("cfi", idx)
  out[["rmsea"]] <- first("rmsea", idx)
  out
}

# The measures psyWiz reports, taken from a raw vector in the chosen variants.
psw_measures <- function(fm, variant) {
  out <- c(chisq = NA_real_, df = NA_real_, pvalue = NA_real_, cfi = NA_real_, tli = NA_real_,
           rmsea = NA_real_, rmsea.lo = NA_real_, rmsea.hi = NA_real_, srmr = NA_real_,
           npar = NA_real_, n = NA_real_)
  if (is.null(fm)) return(out)
  get <- function(nm) if (nm %in% names(fm)) fm[[nm]] else NA_real_
  sfx <- function(v) if (is.na(v)) "" else psw_variant_suffixes[[v]]
  sc <- sfx(variant[["chisq"]]); si <- sfx(variant[["cfi"]]); sr <- sfx(variant[["rmsea"]])
  out[["chisq"]] <- get(paste0("chisq", sc))
  out[["df"]] <- get(paste0("df", sc))
  out[["pvalue"]] <- get(paste0("pvalue", sc))
  out[["cfi"]] <- get(paste0("cfi", si))
  out[["tli"]] <- get(paste0("tli", si))
  out[["rmsea"]] <- get(paste0("rmsea", sr))
  out[["rmsea.lo"]] <- get(paste0("rmsea.ci.lower", sr))
  out[["rmsea.hi"]] <- get(paste0("rmsea.ci.upper", sr))
  out[["srmr"]] <- get("srmr")
  out[["npar"]] <- get("npar")
  out[["n"]] <- get("ntotal")
  out
}

#' Loadings and factor correlations as plain data frames (internal)
#' @keywords internal
#' @noRd
psw_estimates <- function(fit, prep) {
  empty <- list(loadings = NULL, cors = NULL)
  if (!psw_converged(fit)) return(empty)
  pe <- tryCatch(suppressWarnings(lavaan::parameterEstimates(fit, standardized = TRUE, ci = FALSE)),
                 error = function(e) NULL)
  if (is.null(pe)) return(empty)
  if (is.null(pe$group)) pe$group <- 1L
  if (is.null(pe$se)) pe$se <- NA_real_
  if (is.null(pe$pvalue)) pe$pvalue <- NA_real_
  glab <- if (length(prep$groups)) prep$groups[pe$group] else rep(NA_character_, nrow(pe))
  it <- prep$items

  l <- pe$op == "=~"
  loadings <- data.frame(
    group = glab[l],
    factor = unname(prep$fnames[pe$lhs[l]]),
    item = it$name[match(pe$rhs[l], it$id)],
    est = pe$est[l], se = pe$se[l], p = pe$pvalue[l], std = pe$std.all[l],
    stringsAsFactors = FALSE)

  cc <- pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% names(prep$fnames) & pe$rhs %in% names(prep$fnames)
  cors <- data.frame(
    group = glab[cc],
    factor1 = unname(prep$fnames[pe$lhs[cc]]),
    factor2 = unname(prep$fnames[pe$rhs[cc]]),
    r = pe$std.all[cc], p = pe$pvalue[cc],
    stringsAsFactors = FALSE)
  list(loadings = loadings, cors = cors)
}

#' Composite reliability (omega) and AVE per factor from standardised loadings
#'
#' Valid for congeneric factors without cross-loadings or correlated residuals,
#' which is the only kind of model psyWiz builds. With ordinal items the values
#' refer to the latent response variables, not to the observed sum score.
#' @keywords internal
#' @noRd
psw_reliability <- function(loadings) {
  if (is.null(loadings) || nrow(loadings) == 0L) return(NULL)
  key <- paste(loadings$group, loadings$factor, sep = "\r")
  rows <- lapply(split(loadings, factor(key, levels = unique(key))), function(d) {
    s <- sum(d$std)
    data.frame(group = d$group[1L], factor = d$factor[1L],
               omega = s^2 / (s^2 + sum(1 - d$std^2)), ave = mean(d$std^2),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

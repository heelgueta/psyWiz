# engine-syntax.R -- lavaan model syntax for every step of the invariance ladder.
#
# The syntax is written out in full (one modifier per group) instead of relying
# on lavaan's `group.equal`. That keeps the identification constraints explicit,
# lets the ordinal ladder follow Wu & Estabrook (2016), and makes the syntax
# shown to the user exactly what was fitted.
#
# Identification: latent variances fixed to 1 and latent means to 0 in the
# first (reference) group. They are freed in the other groups as soon as the
# loadings (variances) or intercepts (means) are constrained to equality.
#
# Continuous ladder          constrained equal across groups
#   configural               nothing
#   metric                   loadings
#   scalar                   loadings, intercepts
#   strict                   loadings, intercepts, residual variances
#
# Ordinal ladder (Wu & Estabrook, 2016)
#   configural               nothing
#   thresholds               thresholds            (intercepts + scales freed in groups 2+)
#   metric                   thresholds, loadings
#   scalar                   thresholds, loadings, intercepts (fixed to 0)
#   strict                   ... and residual variances (theta parameterisation only)

psw_steps <- function(type, strict = FALSE) {
  s <- if (type == "ordinal") c("configural", "thresholds", "metric", "scalar")
       else c("configural", "metric", "scalar")
  if (strict) c(s, "strict") else s
}

# Each step keeps the constraints of the steps before it.
psw_step_labels <- function(steps, type) {
  lab <- c(configural = "Configural", thresholds = "Thresholds", metric = "Loadings (metric)",
           scalar = "Intercepts (scalar)", strict = "Residuals (strict)")
  unname(lab[steps])
}

# Which parameter type does each step newly constrain?
psw_step_target <- function(step) {
  c(configural = NA, thresholds = "Threshold", metric = "Loading",
    scalar = "Intercept", strict = "Residual variance")[[step]]
}

psw_step_flags <- function(step, type) {
  lvl <- match(step, c("configural", "thresholds", "metric", "scalar", "strict"))
  list(thr = type == "ordinal" && lvl >= 2L, load = lvl >= 3L, int = lvl >= 4L, res = lvl >= 5L)
}

#' Build lavaan syntax for one model (internal)
#'
#' @param prep Output of `psw_prepare()`.
#' @param step One of the steps from `psw_steps()`; ignored for a single group.
#' @param parameterization `"delta"` or `"theta"` (ordinal items only).
#' @return A single string of lavaan syntax using the internal identifiers.
#' @keywords internal
#' @noRd
psw_syntax <- function(prep, step = "configural", parameterization = "delta") {
  it <- prep$items
  G <- length(prep$groups)
  fids <- names(prep$fnames)

  if (G < 2L) {
    out <- vapply(fids, function(f) paste0(f, " =~ ", paste(it$id[it$fid == f], collapse = " + ")), "")
    return(paste(out, collapse = "\n"))
  }

  ord <- prep$type == "ordinal"
  fl <- psw_step_flags(step, prep$type)
  gv <- function(first, rest = first)
    sprintf("c(%s)", paste(c(first, rep(rest, G - 1L)), collapse = ", "))
  out <- character()

  for (f in fids) {
    ids <- it$id[it$fid == f]
    mod <- if (fl$load) vapply(paste0("L_", ids), gv, "") else rep(gv("NA"), length(ids))
    out <- c(out,
             paste0(f, " =~ ", paste0(mod, "*", ids, collapse = " + ")),
             sprintf("%s ~~ %s*%s", f, if (fl$load) gv("1", "NA") else gv("1"), f),
             sprintf("%s ~ %s*1", f, if (fl$int) gv("0", "NA") else gv("0")))
  }

  for (i in it$id) {
    if (!ord) {
      out <- c(out,
               sprintf("%s ~ %s*1", i, if (fl$int) gv(paste0("N_", i)) else gv("NA")),
               sprintf("%s ~~ %s*%s", i, if (fl$res) gv(paste0("E_", i)) else gv("NA"), i))
      next
    }
    k <- seq_len(prep$nthr[[i]])
    tmod <- if (fl$thr) vapply(paste0("T_", i, "_", k), gv, "") else rep(gv("NA"), length(k))
    out <- c(out, paste0(i, " | ", paste0(tmod, "*t", k, collapse = " + ")))
    # intercepts: identified (and freed) once thresholds are equal; fixed to 0 again at the scalar step
    out <- c(out, sprintf("%s ~ %s*1", i, if (fl$thr && !fl$int) gv("0", "NA") else gv("0")))
    # scale / residual variance: freed in groups 2+ once thresholds are equal, but
    # only for items with 3+ categories (not identified for binary items)
    free_scale <- fl$thr && !fl$res && prep$nthr[[i]] >= 2L
    smod <- if (free_scale) gv("1", "NA") else gv("1")
    out <- c(out, if (parameterization == "theta") sprintf("%s ~~ %s*%s", i, smod, i)
                  else sprintf("%s ~*~ %s*%s", i, smod, i))
  }
  paste(out, collapse = "\n")
}

# Syntax with the user's own names, for display only (not parseable if names
# contain spaces or operators).
psw_syntax_display <- function(syntax, prep) {
  map <- c(stats::setNames(prep$items$name, prep$items$id), prep$fnames)
  # translate token by token so a user name that looks like an id is never re-translated
  m <- gregexpr("[A-Za-z_.][A-Za-z0-9_.]*", syntax)
  regmatches(syntax, m) <- lapply(regmatches(syntax, m), function(tok) {
    vapply(tok, function(t) {
      if (t %in% names(map)) return(map[[t]])
      p <- regmatches(t, regexec("^([LNET]_)([if][0-9]+)(_[0-9]+)?$", t))[[1L]]
      if (length(p) && p[3L] %in% names(map)) return(paste0(p[2L], map[[p[3L]]], p[4L]))
      t
    }, "", USE.NAMES = FALSE)
  })
  syntax
}

# engine-robust.R -- how much does an invariance verdict depend on analytic
# choices that are all defensible?
#
# A verdict such as "scalar invariance is not supported" rests on choices the
# analyst could have made differently without being wrong: treating Likert
# items as ordinal or continuous, the estimator, which fit index to look at
# and which cutoff to apply. A multiverse analysis (Steegen et al., 2016)
# makes that dependence visible by running all reasonable combinations.
# invarianceMultiverse() does this for the invariance ladder: it refits the
# ladder under each estimator / item treatment and applies each decision
# rule to every step.

psw_multiverse_rules <- list(
  list(key = "cfi", label = "CFI drop > .010", criteria = "cfi", cutoffs = list(cfi = 0.010)),
  list(key = "cfi002", label = "CFI drop > .002", criteria = "cfi", cutoffs = list(cfi = 0.002)),
  list(key = "rmsea", label = "RMSEA rise > .015", criteria = "rmsea", cutoffs = list(rmsea = 0.015)),
  list(key = "srmr", label = "SRMR rise > .030 / .010", criteria = "srmr", cutoffs = list()),
  list(key = "chisq", label = "Difference test p < .05", criteria = "chisq", cutoffs = list(alpha = 0.05))
)

#' Robustness of invariance verdicts across analytic choices
#'
#' Refits the invariance ladder under every applicable combination of item
#' treatment and estimator (ordinal: WLSMV and ULSMV; continuous: MLR and
#' ML), applies five common decision rules to each step, and reports how
#' often each level of invariance is supported.
#'
#' The decision rules are: CFI drop larger than .010 (Chen, 2007; Cheung &
#' Rensvold, 2002) or than .002 (Meade, Johnson & Braddy, 2008), RMSEA rise
#' larger than .015, SRMR rise larger than .030 (loadings) or .010 (other
#' steps), and a significant chi-square difference test.
#'
#' A step that is supported (or rejected) in nearly all specifications is a
#' robust finding. A step whose verdict flips with the estimator or the rule
#' is one where the data do not settle the question, and should be reported
#' as such.
#'
#' @inheritParams invarianceWiz
#' @param treatments Which item treatments to include: any of `"ordinal"`,
#'   `"continuous"`. Ordinal treatment is skipped (with a note) when the items
#'   cannot be treated as ordinal.
#' @return An object of class `psw_multiverse`: a list with `grid` (one row
#'   per specification, rule and step, with the verdict), `summary` (per step:
#'   number and share of specifications that support it), `ladders` (the
#'   underlying `psw_invariance` objects), `notes` and `info`.
#' @references
#' Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008). Power and
#' sensitivity of alternative fit indices in tests of measurement invariance.
#' *Journal of Applied Psychology, 93*(3), 568-592.
#'
#' Steegen, S., Tuerlinckx, F., Gelman, A., & Vanpaemel, W. (2016).
#' Increasing transparency through a multiverse analysis. *Perspectives on
#' Psychological Science, 11*(5), 702-712.
#' @examples
#' \donttest{
#' data(wizdata)
#' f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
#' invarianceMultiverse(wizdata, f, group = "country")
#' }
#' @export
invarianceMultiverse <- function(data, factors, group, strict = FALSE,
                                 treatments = c("ordinal", "continuous"),
                                 parallel = "auto", ncores = NULL) {
  treatments <- match.arg(treatments, c("ordinal", "continuous"), several.ok = TRUE)
  specs <- list()
  if ("ordinal" %in% treatments)
    specs <- c(specs, list(list(type = "ordinal", estimator = "WLSMV"), list(type = "ordinal", estimator = "ULSMV")))
  if ("continuous" %in% treatments)
    specs <- c(specs, list(list(type = "continuous", estimator = "MLR"), list(type = "continuous", estimator = "ML")))

  notes <- character()
  ladders <- list()
  t0 <- proc.time()[["elapsed"]]
  for (sp in specs) {
    name <- sprintf("%s, %s", if (sp$type == "ordinal") "Ordinal" else "Continuous", sp$estimator)
    lad <- tryCatch(invarianceWiz(data, factors, group, type = sp$type, estimator = sp$estimator,
                                  strict = strict, parallel = parallel, ncores = ncores),
                    error = function(e) conditionMessage(e))
    if (is.character(lad)) {
      notes <- c(notes, sprintf("%s was skipped: %s", name, lad))
      next
    }
    ladders[[name]] <- lad
  }
  if (!length(ladders)) stop("No specification could be fitted. ", paste(notes, collapse = " "), call. = FALSE)

  grid <- do.call(rbind, lapply(names(ladders), function(name) {
    do.call(rbind, lapply(psw_multiverse_rules, function(rule) {
      cp <- invarianceVerdict(ladders[[name]], criteria = rule$criteria, cutoffs = rule$cutoffs)$comparison
      data.frame(specification = name, rule = rule$label, step = cp$step, label = cp$label,
                 verdict = cp$verdict, stringsAsFactors = FALSE)
    }))
  }))
  rownames(grid) <- NULL
  decided <- grid[grid$verdict %in% c("Supported", "Not supported"), , drop = FALSE]
  steps <- unique(grid[c("step", "label")])
  summary <- do.call(rbind, lapply(seq_len(nrow(steps)), function(i) {
    d <- decided[decided$step == steps$step[i], , drop = FALSE]
    data.frame(step = steps$step[i], label = steps$label[i], supported = sum(d$verdict == "Supported"),
               total = nrow(d), share = if (nrow(d)) mean(d$verdict == "Supported") else NA_real_,
               stringsAsFactors = FALSE)
  }))
  summary$reading <- vapply(summary$share, function(s) {
    if (is.na(s)) "not evaluated" else if (s >= 0.9) "robustly supported"
    else if (s <= 0.1) "robustly not supported" else if (s >= 0.5) "mostly supported, depends on choices"
    else "mostly not supported, depends on choices"
  }, "")
  out <- list(grid = grid, summary = summary, ladders = ladders, notes = unique(notes),
              info = list(specifications = names(ladders),
                          rules = vapply(psw_multiverse_rules, function(r) r$label, ""),
                          time = proc.time()[["elapsed"]] - t0))
  class(out) <- "psw_multiverse"
  out
}

# Wide form of the grid: one row per specification and rule, one column per step.
psw_multiverse_wide <- function(x) {
  g <- x$grid
  steps <- unique(g$step)
  key <- paste(g$specification, g$rule, sep = "\r")
  rows <- unique(key)
  out <- g[match(rows, key), c("specification", "rule")]
  for (s in steps) {
    d <- g[g$step == s, , drop = FALSE]
    out[[s]] <- d$verdict[match(rows, paste(d$specification, d$rule, sep = "\r"))]
  }
  rownames(out) <- NULL
  out
}

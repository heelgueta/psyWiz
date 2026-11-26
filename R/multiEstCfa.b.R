# Multi-Estimator CFA. Options, results and panel: jamovi/multiEstCfa.*.yaml;
# jmvtools generates R/multiEstCfa.h.R. The statistics are in cfaWiz()
# (R/engine-multiest.R); this file only moves its output into jamovi tables.
#
# Lifecycle: as in R/jmv-invariance.R. .init() adds one column per ticked
# estimator to the side-by-side tables; .run() fits only when `Run analysis`
# is ticked and caches the result as the state of the `fit` table.

multiEstCfaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "multiEstCfaClass",
  inherit = multiEstCfaBase,
  private = list(
    .estimators = function() {
      opt <- self$options
      psw_all_estimators[vapply(psw_all_estimators, function(e) isTRUE(opt[[paste0("est", e)]]), logical(1))]
    },
    .init = function() {
      est <- private$.estimators()
      res <- self$results
      for (e in est) {
        for (tab in list(res$std, res$unstd, res$cors))
          tab$addColumn(name = e, title = e, type = "number", format = "zto")
      }
      for (e in est)
        res$rel$addColumn(name = paste0("omega_", e), title = e, superTitle = "Composite reliability (\u03C9)",
                          type = "number", format = "zto")
      for (e in est)
        res$rel$addColumn(name = paste0("ave_", e), title = e, superTitle = "AVE",
                          type = "number", format = "zto")
    },
    .run = function() {
      opt <- self$options
      out <- self$results
      factors <- psw_jmv_factors(opt$factors)
      est <- private$.estimators()
      if (!psw_jmv_ready(factors) || !length(est)) {
        psw_jmv_todo(out$todo, c(
          "Drag the items of each factor into its box (at least three items for a single factor) and tick the estimators to compare.",
          "Then tick <b>Run analysis</b>. Nothing is fitted before that, so you can set every option first."))
        return(invisible())
      }

      res <- out$fit$state
      if (is.null(res)) {
        if (!opt$run) {
          prep <- tryCatch(psw_prepare(self$data, factors, type = opt$type),
                           error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
          psw_jmv_todo(out$todo, c(
            sprintf("Ready: %d factor(s), %d items (%s), N = %d. The model will be fitted with %d estimator(s): %s.",
                    length(factors), length(unlist(factors)), prep$type, nrow(prep$data),
                    length(est), paste(est, collapse = ", ")),
            "Tick <b>Run analysis</b> to fit."), prep$notes)
          return(invisible())
        }
        res <- tryCatch(
          cfaWiz(self$data, factors, estimators = est, type = opt$type, missing = opt$missing,
                 parallel = if (opt$parallel) "auto" else FALSE),
          error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
        out$fit$setState(res)
      }
      psw_jmv_todo(out$todo)
      est <- res$info$estimators

      # ---- fit, one row per estimator ----
      f <- res$fit
      psw_jmv_fill(out$fit, lapply(seq_len(nrow(f)), function(i)
        list(estimator = f$estimator[i], treated = f$treated[i], n = f$n[i],
             chisq = f$chisq[i], df = round(f$df[i]), p = f$pvalue[i], cfi = f$cfi[i], tli = f$tli[i],
             rmsea = f$rmsea[i], rmseaLo = f$rmsea.lo[i], rmseaHi = f$rmsea.hi[i], srmr = f$srmr[i])))
      if (isTRUE(opt$highlight)) psw_jmv_flag_fit(out$fit, f)
      for (i in seq_len(nrow(f))) {
        v <- c(chisq = f$v.chisq[i], cfi = f$v.cfi[i], rmsea = f$v.rmsea[i])
        if (!anyNA(v) && any(v != "standard"))
          out$fit$addFootnote(rowNo = i, col = "estimator", psw_jmv_variant_text(v))
      }
      out$fit$setNote("info", paste(
        "ML, MLR and MLM treat the items as continuous; the least-squares estimators treat ordinal items as ordinal (polychoric correlations).",
        psw_jmv_time_text(res$info)), init = FALSE)
      psw_jmv_notes(out$fit, res$notes)

      # ---- loadings, one column per estimator ----
      loadings <- function(table, value) {
        w <- psw_wide(res$loadings, c("factor", "item"), value, est)
        if (is.null(w)) return()
        if (opt$sortLoadings) {
          s <- psw_wide(res$loadings, c("factor", "item"), "std", est)
          avg <- rowMeans(abs(as.matrix(s[est])), na.rm = TRUE)
          w <- w[order(match(w$factor, names(res$info$factors)), -avg), , drop = FALSE]
        }
        psw_jmv_fill(table, lapply(seq_len(nrow(w)), function(i) as.list(w[i, , drop = FALSE])))
        if (opt$highlight && value == "std") {
          for (i in seq_len(nrow(w))) for (e in est)
            if (isTRUE(abs(w[[e]][i]) < 0.40))
              table$addFormat(rowNo = i, col = e, format = jmvcore::Cell.NEGATIVE)
          table$setNote("hl", "Loadings below 0.40 (in absolute value) are shown in red.", init = FALSE)
        }
      }
      if (opt$showStd) loadings(out$std, "std")
      if (opt$showUnstd) loadings(out$unstd, "est")

      if (opt$showCors) {
        w <- psw_wide(res$cors, c("factor1", "factor2"), "r", est)
        if (!is.null(w))
          psw_jmv_fill(out$cors, lapply(seq_len(nrow(w)), function(i) as.list(w[i, , drop = FALSE])))
        else out$cors$setNote("one", "The model has a single factor.", init = FALSE)
      }

      if (opt$showRel && !is.null(res$reliability)) {
        om <- psw_wide(res$reliability, "factor", "omega", est)
        av <- psw_wide(res$reliability, "factor", "ave", est)
        psw_jmv_fill(out$rel, lapply(seq_len(nrow(om)), function(i)
          c(list(factor = om$factor[i]),
            stats::setNames(as.list(om[i, est, drop = FALSE]), paste0("omega_", est)),
            stats::setNames(as.list(av[i, est, drop = FALSE]), paste0("ave_", est)))))
        if (isTRUE(opt$highlight)) {
          for (i in seq_len(nrow(om))) for (e in est) {
            if (isTRUE(om[[e]][i] < 0.70))
              out$rel$addFormat(rowNo = i, col = paste0("omega_", e), format = jmvcore::Cell.NEGATIVE)
            if (isTRUE(av[[e]][i] < 0.50))
              out$rel$addFormat(rowNo = i, col = paste0("ave_", e), format = jmvcore::Cell.NEGATIVE)
          }
        }
        out$rel$setNote("what", paste0(
          "Computed from the standardized loadings. For estimators that treat the items as ordinal, ",
          "the values refer to the latent response variables, not to the observed sum score.",
          if (isTRUE(opt$highlight)) " In red: \u03C9 below .70 or AVE below .50." else ""), init = FALSE)
      }

      invisible()
    }
  )
)

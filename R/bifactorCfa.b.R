# Bifactor CFA. Options, results and panel: jamovi/bifactorCfa.*.yaml (generated
# by tools/gen_yaml.py); jmvtools generates R/bifactorCfa.h.R. The statistics
# are in bifactorWiz() (R/engine-bifactor.R). Lifecycle as in
# R/jmv-invariance.R: nothing is fitted until `Run analysis` is ticked, and
# the result is cached as the state of the `fit` table.

bifactorCfaClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "bifactorCfaClass",
  inherit = bifactorCfaBase,
  private = list(
    .run = function() {
      opt <- self$options
      out <- self$results
      factors <- psw_jmv_factors(opt$factors)
      if (length(factors) < 2L || any(lengths(factors) < 3L)) {
        psw_jmv_todo(out$todo, c(
          "Drag the items of each <b>specific</b> factor (subscale) into its box: at least two factors with three or more items each. Every item will also load on the general factor.",
          "Then tick <b>Run analysis</b>. Nothing is fitted before that, so you can set every option first."))
        return(invisible())
      }
      general <- trimws(opt$generalName)
      if (!nzchar(general)) general <- "General"
      res <- out$fit$state
      if (is.null(res)) {
        if (!opt$run) {
          prep <- tryCatch(psw_prepare(self$data, factors, type = opt$type),
                           error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
          psw_jmv_todo(out$todo, c(
            sprintf("Ready: one general and %d specific factors, %d items treated as <b>%s</b>, N = %d.",
                    length(factors), length(unlist(factors)), prep$type, nrow(prep$data)),
            "Tick <b>Run analysis</b> to fit the model."), prep$notes)
          return(invisible())
        }
        res <- tryCatch(
          bifactorWiz(self$data, factors, general = general, type = opt$type, estimator = opt$estimator,
                      missing = opt$missing),
          error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
        out$fit$setState(res)
      }
      psw_jmv_todo(out$todo)
      info <- res$info

      psw_jmv_fill(out$fit, psw_jmv_fit_rows(res$fit))
      if (isTRUE(opt$highlight)) psw_jmv_flag_fit(out$fit, res$fit)
      out$fit$setNote("info", sprintf(
        "Bifactor: every item loads on %s and on its specific factor; all factors uncorrelated. N = %d. Estimator: %s; items treated as %s; missing values: %s. %s lavaan %s.",
        info$general, info$n, info$estimator, info$type, info$missing, paste0(psw_jmv_variant_text(info$variant), ""),
        info$lavaan), init = FALSE)
      out$fit$setNote("fit", paste0(
        "A bifactor model has more parameters and tends to fit better than the correlated-factors model even ",
        "when it is not the better description; read the indices below, not only the fit."), init = FALSE)
      psw_jmv_notes(out$fit, res$notes)

      if (opt$showIndices) {
        ix <- res$indices
        psw_jmv_fill(out$indices, lapply(seq_len(nrow(ix)), function(i) as.list(ix[i, , drop = FALSE])))
        out$indices$setNote("what", sprintf(paste0(
          "\u03C9: reliability of the score. \u03C9 hierarchical: for the total score, the share of its variance due to %s; ",
          "for a subscale, the share due to its specific factor once %s is removed (low values mean the subscale ",
          "score says little beyond the total). ECV: share of the explained variance that is general. ",
          "%.0f%% of the item pairs come from different specific factors (PUC)."),
          info$general, info$general, 100 * info$puc), init = FALSE)
      }
      if (opt$showLoadings) {
        l <- res$loadings
        psw_jmv_fill(out$loadings, lapply(seq_len(nrow(l)), function(i) as.list(l[i, , drop = FALSE])))
        if (isTRUE(opt$highlight)) {
          for (i in seq_len(nrow(l))) {
            if (isTRUE(l$unique[i] > 0.80 || l$unique[i] < 0))
              out$loadings$addFormat(rowNo = i, col = "unique", format = jmvcore::Cell.NEGATIVE)
            if (isTRUE(abs(l$general[i]) < 0.30 && abs(l$own[i]) < 0.30)) {
              out$loadings$addFormat(rowNo = i, col = "general", format = jmvcore::Cell.NEGATIVE)
              out$loadings$addFormat(rowNo = i, col = "own", format = jmvcore::Cell.NEGATIVE)
            }
          }
        }
        out$loadings$setNote("what", paste0(
          "Standardized loadings. Unique: variance explained by neither factor. Item ECV: share of the item's ",
          "explained variance that is general.",
          if (isTRUE(opt$highlight)) " In red: high uniqueness (> .80) or negative uniqueness, and items weak (< .30) on both factors." else ""),
          init = FALSE)
      }

      fit_text <- if (!is.null(opt$plotFit) && opt$plotFit != "none" && !is.null(res$fit)) {
        psw_format_fit_line(res$fit[1L, , drop = FALSE], info$estimator)
      } else NULL
      fit_pos <- if (!is.null(opt$plotFit)) opt$plotFit else "none"
      show_guide <- !isFALSE(opt$showPlotGuide)

      if (opt$showAntenna) {
        psw_jmv_plot(out$antenna, "antenna", list(L = res$loadings, general = info$general,
                                                  fit_text = fit_text, fit_pos = fit_pos))
        if (!is.null(out$antennaDesc) && show_guide) out$antennaDesc$setContent(psw_jmv_antenna_desc())
      }
      if (opt$showFacing) {
        psw_jmv_plot(out$facing, "facing", list(L = res$loadings, general = info$general,
                                                fit_text = fit_text, fit_pos = fit_pos))
        if (!is.null(out$facingDesc) && show_guide) out$facingDesc$setContent(psw_jmv_facing_desc())
      }
      if (opt$showPath)
        psw_jmv_plot(out$pathFig, "bfpath", list(L = res$loadings, general = info$general,
                                                 numbers = opt$pathNumbers, width = opt$pathWidth,
                                                 fit_text = fit_text, fit_pos = fit_pos))
      invisible()
    },
    .diagram = function(image, ...) psw_jmv_diagram(image)
  )
)

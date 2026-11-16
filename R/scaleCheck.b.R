# Scale Check. Options, results and panel: jamovi/scaleCheck.*.yaml (generated
# by tools/gen_yaml.py); jmvtools generates R/scaleCheck.h.R. The statistics
# are in scaleWiz() (R/engine-scale.R). Lifecycle as in R/jmv-invariance.R:
# nothing is fitted until `Run analysis` is ticked, and the result is cached
# as the state of the `reliability` table (the one table that is always
# shown). By default only reliability and the item table are shown; the
# tau-equivalence ladder and the Guttman / Mokken statistics are computed
# anyway (they are cheap for one scale) and shown on request.

scaleCheckClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "scaleCheckClass",
  inherit = scaleCheckBase,
  private = list(
    .run = function() {
      opt <- self$options
      out <- self$results
      vars <- opt$vars
      if (length(vars) < 3L) {
        psw_jmv_todo(out$todo, c(
          "Choose the items of one scale (at least three).",
          "Then tick <b>Run analysis</b>. Nothing is fitted before that, so you can set every option first."))
        return(invisible())
      }
      name <- trimws(opt$scaleName)
      if (!nzchar(name)) name <- "Scale"

      res <- out$reliability$state
      if (is.null(res)) {
        if (!opt$run) {
          prep <- tryCatch(psw_prepare(self$data, stats::setNames(list(vars), name), type = opt$type),
                           error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
          psw_jmv_todo(out$todo, c(
            sprintf("Ready: %d items treated as <b>%s</b>, N = %d.", length(vars), prep$type, nrow(prep$data)),
            "Tick <b>Run analysis</b> to compute."), prep$notes)
          return(invisible())
        }
        res <- tryCatch(
          scaleWiz(self$data, vars, name = name, type = opt$type, estimator = opt$estimator,
                   missing = opt$missing, parallel = if (opt$parallel) "auto" else FALSE),
          error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
        out$reliability$setState(res)
      }
      psw_jmv_todo(out$todo)
      cr <- psw_jmv_criteria(opt)
      res <- invarianceVerdict(res, criteria = cr$criteria, cutoffs = cr$cutoffs)
      info <- res$info

      # ---- reliability ----
      rl <- res$reliability
      psw_jmv_fill(out$reliability, lapply(seq_len(nrow(rl)), function(i)
        list(coefficient = rl$coefficient[i], value = rl$value[i], basis = rl$basis[i])))
      out$reliability$setNote("what", paste0(
        "Alpha equals reliability only if the items are tau-equivalent (see More); lambda-2 is never lower than ",
        "alpha and is still a lower bound; omega comes from the congeneric model and does not need equal loadings. ",
        if (info$type == "ordinal")
          "The ordinal coefficients describe the latent response variables, not the sum score people actually compute. "
        else "",
        sprintf("N = %d complete cases. Estimator for omega: %s.", info$n.complete, info$estimator)), init = FALSE)
      psw_jmv_notes(out$reliability, res$notes)

      if (opt$tau) {
        # ---- models and their comparison ----
        psw_jmv_fill(out$fit, psw_jmv_fit_rows(res$fit))
        out$fit$setNote("info", sprintf(
          "Congeneric: every item has its own loading. Tau-equivalent: equal loadings. %sEstimator: %s; items treated as %s; missing values: %s. %s %s",
          if (info$type == "ordinal") "" else "Parallel: equal loadings and equal residual variances. ",
          info$estimator, info$type, info$missing, paste0(psw_jmv_variant_text(info$variant), ""),
          psw_jmv_time_text(info)), init = FALSE)

        cp <- res$comparison
        psw_jmv_fill(out$comparison, lapply(seq_len(nrow(cp)), function(i)
          list(model = cp$label[i], vs = cp$vs[i], dchisq = cp$dchisq[i], ddf = cp$ddf[i], p = cp$p[i],
               dcfi = cp$dcfi[i], drmsea = cp$drmsea[i], dsrmr = cp$dsrmr[i], verdict = cp$verdict[i])))
        flagcols <- c(cfi = "dcfi", rmsea = "drmsea", srmr = "dsrmr", chisq = "p")
        if (opt$highlight) for (i in seq_len(nrow(cp))) for (k in cr$criteria)
          if (isTRUE(cp[[paste0("flag.", k)]][i]))
            out$comparison$addFormat(rowNo = i, col = flagcols[[k]], format = jmvcore::Cell.NEGATIVE)
        reading <- psw_scale_reading(res)
        for (i in seq_along(reading)) out$comparison$setNote(paste0("read", i), reading[i], init = FALSE)
        if (info$type == "ordinal")
          out$comparison$setNote("ord", paste0(
            "With ordinal items the latent responses have unit variance, so equal loadings imply equal residual ",
            "variances: tau-equivalence and parallelism are the same model."), init = FALSE)
      }

      # ---- items ----
      if (opt$showItems) {
        it <- res$items
        psw_jmv_fill(out$itemStats, lapply(seq_len(nrow(it)), function(i)
          list(item = it$item[i], mean = it$mean[i], sd = it$sd[i], loading = it$loading[i], r2 = it$r2[i],
               restR = it$rest.r[i], alphaDrop = it$alpha.drop[i], H = it$H[i])))
        alpha <- rl$value[rl$coefficient == "Cronbach's alpha"]
        if (opt$highlight) for (i in seq_len(nrow(it))) {
          if (isTRUE(abs(it$loading[i]) < 0.40)) out$itemStats$addFormat(rowNo = i, col = "loading", format = jmvcore::Cell.NEGATIVE)
          if (isTRUE(it$H[i] < 0.30)) out$itemStats$addFormat(rowNo = i, col = "H", format = jmvcore::Cell.NEGATIVE)
          if (isTRUE(it$alpha.drop[i] > alpha)) out$itemStats$addFormat(rowNo = i, col = "alphaDrop", format = jmvcore::Cell.NEGATIVE)
        }
        out$itemStats$setNote("what", paste0(
          "Loading and R\u00B2: standardized, from the one-factor model. Item-rest r: correlation with the sum of the ",
          "other items.",
          if (opt$showScal) " H: item scalability (Loevinger); below .30 the item does not order people like the others do." else "",
          if (opt$highlight) " In red: loading below .40, alpha higher without the item, or H below .30." else ""),
          init = FALSE)
      }

      # ---- scalability ----
      if (opt$showScal) {
        sc <- res$guttman[[opt$guttmanOrder]]$table
        psw_jmv_fill(out$scalability, lapply(seq_len(nrow(sc)), function(i)
          list(statistic = sc$statistic[i], value = sc$value[i], reading = sc$reading[i])))
        out$scalability$setNote("what", paste0(
          "A Guttman scale is one where anyone who endorses a hard item also endorses all easier ones. ",
          "Loevinger's H measures how close the data come to that (1 = perfect, 0 = unrelated items) and works for ",
          "ordinal items too, whatever the order of the items. CR, MMR and CS are the classical scalogram statistics; ",
          "they depend on the order the items are tested in. ",
          if (opt$guttmanOrder == "entered")
            "Here that is the order you entered them (easiest first), which tests your ordering; for items with more than two categories the item steps are sorted by endorsement instead."
          else "Here the items are sorted from most to least endorsed, the order that fits this sample best."),
          init = FALSE)
      }

      if (opt$showScalogram) {
        psw_jmv_plot(out$scalogram, "scalogram", list(sc = res$guttman[[opt$guttmanOrder]]$scalogram))
      }

      if (!is.null(res$loadings)) {
        l <- res$loadings
        l$group <- sprintf("N = %d", info$n)
        psw_jmv_figures(out, opt, l, NULL)
      }
      invisible()
    },
    .diagram = function(image, ...) psw_jmv_diagram(image)
  )
)

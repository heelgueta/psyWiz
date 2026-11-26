# jmv-invariance.R -- jamovi implementation of Group CFA (groupCfa). Up to 0.3
# this was shared by two analyses (one-factor and correlated-factors), which
# is why it takes `factors` as an argument instead of reading the options.
#
# Lifecycle
#   .init()  adds one loadings column per group (the groups are known from
#            the data header) and decides which tables are visible.
#   .run()   1. not enough specified  -> instructions
#            2. `Run analysis` unticked and nothing cached -> group sizes and
#               the plan, no fitting
#            3. otherwise fit (invarianceWiz) or reuse the cached result. The
#               cache is the state of the `fit` table, cleared by jamovi when
#               an option in its clearWith changes.
#            4. apply the decision criteria (cheap, never refits) and fill
#               the tables.
#            5. on request, the partial-invariance search and the robustness
#               check, each with its own cache (the state of its first table)
#               because each refits models.

psw_jmv_inv_init <- function(self, factors) {
  opt <- self$options
  self$results$culprits$setVisible(opt$culprits != "none")
  g <- opt$group
  if (!is.null(g) && g %in% names(self$data)) {
    lv <- levels(droplevels(as.factor(self$data[[g]])))
    if (!length(lv)) lv <- levels(as.factor(self$data[[g]]))
    for (k in seq_along(lv))
      self$results$loadings$addColumn(name = paste0("g", k), title = lv[k], type = "number", format = "zto")
  }
}

psw_jmv_inv_run <- function(self, factors) {
  opt <- self$options
  out <- self$results
  if (!psw_jmv_ready(factors) || is.null(opt$group)) {
    psw_jmv_todo(out$todo, c(
      "Drag the items of each factor into its box (at least three items for a single factor) and choose a grouping variable.",
      "Then tick <b>Run analysis</b>. Nothing is fitted before that, so you can set every option first."))
    return(invisible())
  }

  res <- out$fit$state
  if (is.null(res)) {
    prep <- tryCatch(psw_prepare(self$data, factors, group = opt$group, type = opt$type),
                     error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
    psw_jmv_fill(out$groups, lapply(seq_along(prep$n), function(i)
      list(group = names(prep$n)[i], n = as.integer(prep$n[[i]]))))
    if (!opt$run) {
      steps <- psw_step_labels(psw_steps(prep$type, opt$strict), prep$type)
      est <- if (opt$estimator == "auto") (if (prep$type == "ordinal") "WLSMV" else "MLR") else opt$estimator
      psw_jmv_todo(out$todo, c(
        sprintf("Ready: %d factor(s), %d items treated as <b>%s</b>, %d groups. %d models will be fitted with %s: %s.",
                length(factors), length(unlist(factors)), prep$type, length(prep$groups),
                length(steps), est, paste(steps, collapse = ", ")),
        "Tick <b>Run analysis</b> to fit them."), prep$notes)
      return(invisible())
    }
    res <- tryCatch(
      invarianceWiz(self$data, factors, group = opt$group, type = opt$type,
                    estimator = opt$estimator, missing = opt$missing, strict = opt$strict,
                    parallel = if (opt$parallel) "auto" else FALSE),
      error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
    # the prepared data would be saved with the state (and in the .omv file); it is rebuilt when needed
    res$prep <- NULL
    out$fit$setState(res)
  }
  psw_jmv_todo(out$todo)

  criteria <- c("cfi", "rmsea", "srmr", "chisq")[c(opt$useCfi, opt$useRmsea, opt$useSrmr, opt$useChisq)]
  res <- invarianceVerdict(res, criteria = criteria,
                           cutoffs = list(cfi = opt$cutCfi, rmsea = opt$cutRmsea,
                                          srmr.loadings = opt$cutSrmrLoad, srmr = opt$cutSrmr,
                                          alpha = opt$alpha))

  # ---- groups ----
  psw_jmv_fill(out$groups, lapply(seq_len(nrow(res$groups)), function(i)
    list(group = res$groups$group[i], n = res$groups$n[i])))

  # ---- model fit ----
  f <- res$fit
  psw_jmv_fill(out$fit, lapply(seq_len(nrow(f)), function(i)
    list(model = f$label[i], chisq = f$chisq[i], df = round(f$df[i]), p = f$pvalue[i],
         cfi = f$cfi[i], tli = f$tli[i], rmsea = f$rmsea[i],
         rmseaLo = f$rmsea.lo[i], rmseaHi = f$rmsea.hi[i], srmr = f$srmr[i])))
  if (isTRUE(opt$highlight)) psw_jmv_flag_fit(out$fit, f)
  info <- res$info
  out$fit$setNote("info", sprintf(
    "Estimator: %s. Items treated as %s. Missing values: %s. Each model keeps the constraints of the models above it. %s %s",
    info$estimator, info$type, info$missing,
    paste0(psw_jmv_variant_text(info$variant), ""), psw_jmv_time_text(info)), init = FALSE)
  psw_jmv_notes(out$fit, res$notes)

  # ---- step comparisons ----
  cp <- res$comparison
  psw_jmv_fill(out$comparison, lapply(seq_len(nrow(cp)), function(i)
    list(model = cp$label[i], vs = cp$vs[i], dchisq = cp$dchisq[i], ddf = cp$ddf[i], p = cp$p[i],
         dcfi = cp$dcfi[i], drmsea = cp$drmsea[i], dsrmr = cp$dsrmr[i], verdict = cp$verdict[i])))
  flagcols <- c(cfi = "dcfi", rmsea = "drmsea", srmr = "dsrmr", chisq = "p")
  if (opt$highlight) for (i in seq_len(nrow(cp))) for (k in criteria)
    if (isTRUE(cp[[paste0("flag.", k)]][i]))
      out$comparison$addFormat(rowNo = i, col = flagcols[[k]], format = jmvcore::Cell.NEGATIVE)
  cut <- res$cutoffs
  rule <- c(cfi = sprintf("CFI drops by more than %s", format(cut$cfi)),
            rmsea = sprintf("RMSEA rises by more than %s", format(cut$rmsea)),
            srmr = sprintf("SRMR rises by more than %s (loadings step) or %s (other steps)",
                           format(cut$srmr.loadings), format(cut$srmr)),
            chisq = sprintf("the \u03C7\u00B2 difference test has p < %s", format(cut$alpha)))
  out$comparison$setNote("rule", if (length(criteria)) paste0(
    "Each step is compared with the previous one. A value is flagged",
    if (opt$highlight) " (shown in red)" else "", " when ",
    paste(rule[criteria], collapse = "; "),
    ". Invariance is 'Supported' when nothing is flagged, 'Not supported' when every selected criterion is flagged, ",
    "and 'Ambiguous' when the criteria disagree. Default cutoffs follow Chen (2007).")
    else "No decision criterion is selected (see Invariance Cutoffs).", init = FALSE)
  if (info$type == "ordinal")
    out$comparison$setNote("ladder", paste0(
      "Ordinal items: steps and identification constraints follow Wu & Estabrook (2016). ",
      "For items with three categories the Thresholds model is equivalent to the Configural model."), init = FALSE)

  # ---- worst-fitting constraints ----
  if (opt$culprits != "none") {
    steps <- if (opt$culprits == "all") cp$step
             else cp$step[cp$verdict %in% c("Not supported", "Ambiguous")]
    cu <- res$culprits
    rows <- list()
    for (s in steps) {
      d <- utils::head(cu[cu$step == s, , drop = FALSE], opt$culpritsN)
      for (i in seq_len(nrow(d)))
        rows[[length(rows) + 1L]] <- list(
          model = d$label[i], parameter = d$parameter[i],
          item = if (nzchar(d$detail[i])) sprintf("%s, threshold %s", d$item[i], sub("^t", "", d$detail[i])) else d$item[i],
          group = d$group[i], x2 = d$x2[i], df = d$df[i], p = d$p[i])
    }
    psw_jmv_fill(out$culprits, rows)
    out$culprits$setNote("what", paste0(
      if (!length(rows)) (if (opt$culprits == "failed") "No step was flagged, so there is nothing to list. " else "No score tests are available. ") else "",
      "Score tests: the expected drop in the model \u03C7\u00B2 if this single constraint were released for the group shown ",
      "(relative to the first group). They are not robust-corrected; use them to rank constraints, not as exact tests."),
      init = FALSE)
  }

  # ---- partial invariance: its own cache, on the partial ladder table ----
  # jamovi clears that state when the model or a decision criterion changes,
  # because the search depends on the verdicts.
  if (opt$partial) {
    part <- out$partialLadder$state
    if (is.null(part) && opt$run) {
      part <- tryCatch({
        res$prep <- psw_prepare(self$data, factors, group = opt$group, type = info$type)
        invariancePartial(res, max.free = opt$maxFree)
      }, error = function(e) conditionMessage(e))
      res$prep <- NULL
      if (is.list(part)) out$partialLadder$setState(part)
    }
    if (is.null(part)) {
      out$partialLadder$setNote("todo", "Tick 'Run analysis' to run the search.", init = FALSE)
    } else if (is.character(part)) {
      out$partialLadder$setNote("err", paste("Partial invariance could not be computed:", part), init = FALSE)
    } else {
      pcp <- part$comparison
      psw_jmv_fill(out$partialLadder, lapply(seq_len(nrow(pcp)), function(i)
        list(model = pcp$label[i], vs = pcp$vs[i], dchisq = pcp$dchisq[i], ddf = pcp$ddf[i], p = pcp$p[i],
             dcfi = pcp$dcfi[i], drmsea = pcp$drmsea[i], dsrmr = pcp$dsrmr[i],
             verdict = pcp$verdict[i], freed = pcp$freed[i])))
      if (opt$highlight) for (i in seq_len(nrow(pcp))) for (k in criteria)
        if (isTRUE(pcp[[paste0("flag.", k)]][i]))
          out$partialLadder$addFormat(rowNo = i, col = flagcols[[k]], format = jmvcore::Cell.NEGATIVE)
      out$partialLadder$setNote("summary", part$summary, init = FALSE)
      out$partialLadder$setNote("how", paste0(
        "The ladder refitted: at a step that is not supported, the item with the largest score test is left free ",
        "in every group and the step is judged again, keeping at least two constrained items per factor. ",
        "Thresholds are never released. This is a search guided by the data: it describes where the misfit ",
        "sits in this sample."), init = FALSE)
      psw_jmv_notes(out$partialLadder, part$notes)
      rel <- part$released
      if (is.null(rel))
        out$partialReleased$setNote("none", "Nothing was released: every step was supported, or no release was allowed.", init = FALSE)
      else
        psw_jmv_fill(out$partialReleased, lapply(seq_len(nrow(rel)), function(i)
          list(model = rel$step[i], parameter = rel$parameter[i], item = rel$item[i], factor = rel$factor[i],
               group = rel$group[i], x2 = rel$x2[i], df = rel$df[i], p = rel$p[i], reason = rel$reason[i])))
    }
  }

  # ---- loadings of the configural model, one column per group ----
  if (opt$showLoadings && !is.null(res$loadings)) {
    l <- res$loadings[res$loadings$step == "configural", , drop = FALSE]
    if (nrow(l)) {
      key <- paste(l$factor, l$item, sep = "\r")
      rows <- lapply(unique(key), function(k) {
        d <- l[key == k, , drop = FALSE]
        vals <- as.list(stats::setNames(d$std[match(res$groups$group, d$group)],
                                        paste0("g", seq_len(nrow(res$groups)))))
        c(list(factor = d$factor[1L], item = d$item[1L]), vals)
      })
      psw_jmv_fill(out$loadings, rows)
      if (isTRUE(opt$highlight)) {
        gcols <- paste0("g", seq_len(nrow(res$groups)))
        for (i in seq_along(rows)) for (gc in gcols)
          if (isTRUE(abs(rows[[i]][[gc]]) < 0.40))
            out$loadings$addFormat(rowNo = i, col = gc, format = jmvcore::Cell.NEGATIVE)
        out$loadings$setNote("hl", "In red: standardized loadings below .40.", init = FALSE)
      }
    }
  }

  # ---- latent mean comparisons ----
  if (opt$latentMeans || opt$latentMeansPlot) {
    ref_grp <- if (nzchar(opt$latentMeansRef)) trimws(opt$latentMeansRef) else NULL
    if (isTRUE(opt$partial) && !is.null(out$partialLadder$state) && is.list(out$partialLadder$state)) {
      res$partial <- out$partialLadder$state
    }
    means_obj <- tryCatch(
      invarianceMeans(res, reference = ref_grp),
      error = function(e) conditionMessage(e)
    )

    if (is.character(means_obj)) {
      out$latentMeans$setNote("err", paste("Latent means could not be computed:", means_obj), init = FALSE)
    } else if (isTRUE(opt$latentMeans)) {
      if (!isTRUE(attr(means_obj, "supported"))) {
        out$latentMeans$setNote("unsupported",
          "Latent means cannot be compared: scalar invariance is not supported across groups. Scalar invariance requires equal loadings and intercepts (or thresholds), ensuring that group differences in factor means reflect true construct differences rather than item bias.",
          init = FALSE)
      } else {
        psw_jmv_fill(out$latentMeans, lapply(seq_len(nrow(means_obj)), function(i) {
          list(
            factor = means_obj$factor[i],
            group = means_obj$group[i],
            mean = means_obj$mean[i],
            se = means_obj$se[i],
            ci_lower = means_obj$ci_lower[i],
            ci_upper = means_obj$ci_upper[i],
            z = means_obj$z[i],
            p = means_obj$p[i],
            d_ref = means_obj$d_ref[i],
            d_pool = means_obj$d_pool[i]
          )
        }))
        ref_text <- sprintf("Differences relative to reference group '%s' (mean = 0.00).", attr(means_obj, "reference"))
        type_text <- if (identical(attr(means_obj, "model_type"), "partial"))
          sprintf("Estimated under partial scalar invariance (freed items: %s).", paste(attr(means_obj, "freed_items"), collapse = ", "))
        else "Estimated under full scalar invariance."
        eff_text <- "Cohen's d (ref) standardizes by the reference group factor SD; Cohen's d (pool) standardizes by the pooled factor SD (Hancock, 2001)."
        out$latentMeans$setNote("summary", paste(ref_text, type_text, eff_text), init = FALSE)
      }
    }
    if (isTRUE(opt$latentMeansPlot) && !is.character(means_obj)) {
      psw_jmv_plot(out$latentMeansPlot, "means", list(means = means_obj))
    }
  }

  # ---- figure: the render function (.diagram in the *.b.R files) draws from this state ----
  if ((opt$showComb || opt$showPath) && !is.null(res$loadings)) {
    l <- res$loadings[res$loadings$step == "configural", , drop = FALSE]
    cc <- if (is.null(res$cors)) NULL else res$cors[res$cors$step == "configural", , drop = FALSE]
    f_conf <- if (!is.null(res$fit) && any(res$fit$step == "configural")) res$fit[res$fit$step == "configural", , drop = FALSE] else NULL
    if (nrow(l)) psw_jmv_figures(out, opt, l, cc, caption = "Configural model.", fit = f_conf, estimator = info$estimator)
  }

  # ---- robustness of the verdict: its own cache, on the summary table ----
  if (opt$robust) {
    mv <- out$robustSummary$state
    if (is.null(mv) && opt$run) {
      mv <- tryCatch(invarianceMultiverse(self$data, factors, group = opt$group, strict = opt$strict,
                                          parallel = if (opt$parallel) "auto" else FALSE),
                     error = function(e) conditionMessage(e))
      if (is.list(mv)) {
        mv$ladders <- NULL
        out$robustSummary$setState(mv)
      }
    }
    if (is.null(mv)) {
      out$robustSummary$setNote("todo", "Tick 'Run analysis' to compute the robustness check.", init = FALSE)
    } else if (is.character(mv)) {
      out$robustSummary$setNote("err", paste("The robustness check could not be computed:", mv), init = FALSE)
    } else {
      sm <- mv$summary
      psw_jmv_fill(out$robustSummary, lapply(seq_len(nrow(sm)), function(i)
        list(model = sm$label[i], supported = sprintf("%d of %d", sm$supported[i], sm$total[i]),
             share = sm$share[i], reading = sm$reading[i])))
      out$robustSummary$setNote("what", paste0(
        "The ladder was refitted under each combination of item treatment and estimator (",
        paste(mv$info$specifications, collapse = "; "), ") and five decision rules were applied to every step. ",
        "A verdict that holds in nearly all specifications is robust; one that flips depends on choices ",
        "the data cannot settle, and should be reported as uncertain."), init = FALSE)
      psw_jmv_notes(out$robustSummary, mv$notes)
      w <- psw_multiverse_wide(mv)
      psw_jmv_fill(out$robustGrid, lapply(seq_len(nrow(w)), function(i) {
        row <- list(specification = w$specification[i], rule = w$rule[i])
        for (s in c("thresholds", "metric", "scalar", "strict"))
          row[[s]] <- if (s %in% names(w)) w[[s]][i] else NA
        row
      }))
      out$robustGrid$setNote("thr", "The thresholds step exists only when items are treated as ordinal.", init = FALSE)
    }
  }
  invisible()
}


# jmv-invariance.R -- jamovi implementation shared by One-Factor Group CFA
# (groupCfaOne) and Correlated-Factors Group CFA (groupCfaCorr). Their *.b.R
# files only translate their own model options into a `factors` list.
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
    else "No decision criterion is selected (see Decision Criteria).", init = FALSE)
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
    }
  }

  # ---- diagram: the render function (.diagram in the *.b.R files) draws from this state ----
  if (opt$showDiagram && !is.null(res$loadings)) {
    l <- res$loadings[res$loadings$step == "configural", , drop = FALSE]
    cc <- if (is.null(res$cors)) NULL else res$cors[res$cors$step == "configural", , drop = FALSE]
    if (nrow(l)) {
      size <- psw_diagram_size(l, cc, panel_col = "group")
      out$diagram$setSize(size[["width"]], size[["height"]])
      out$diagram$setState(list(loadings = l, cors = cc, highlight = opt$highlight))
    }
  }
  invisible()
}

# Render function shared by both group analyses.
psw_jmv_diagram <- function(image) {
  st <- image$state
  if (is.null(st)) return(FALSE)
  print(psw_diagram(st$loadings, st$cors, panel_col = "group", highlight = isTRUE(st$highlight)))
  TRUE
}

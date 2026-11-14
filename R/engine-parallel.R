# engine-parallel.R -- run several independent model fits, in parallel when it pays.
#
# A single lavaan fit is single-threaded and cannot be split up without
# rewriting lavaan's optimiser. What *can* run concurrently are the separate
# models an analysis needs (the steps of an invariance ladder, or one model per
# estimator). Strategy for parallel = "auto":
#
#   1. Fit the models one after the other in this process, timing each.
#   2. As soon as one of them was slow and at least two models remain, start
#      background R workers (PSOCK, so it works on Windows, macOS and Linux)
#      and fit the rest concurrently. If workers cannot be started we
#      silently carry on as in step 1.
#
# Starting workers costs a second or two, which is why quick models stay
# sequential.
#
# Warm starts (starting each invariance step from the previous step's
# estimates) were tried and dropped: they saved no measurable time, because
# the cost is dominated by the sample statistics, and lavaan::modindices()
# returned zeros for models fitted with `start = <fitted model>`.

psw_run_task <- function(task) {
  switch(task$kind,
         inv = psw_inv_task(task),
         est = psw_est_task(task),
         stop("unknown task kind"))
}

psw_default_cores <- function(ntasks) {
  n <- suppressWarnings(parallel::detectCores())
  if (is.na(n)) n <- 1L
  n <- n - 1L
  # CRAN's check machines allow at most two cores
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) n <- min(n, 2L)
  max(1L, min(ntasks, n))
}

# Start workers that can see the same packages as this process; NULL on failure.
psw_cluster <- function(n) {
  cl <- tryCatch(suppressWarnings(parallel::makePSOCKcluster(n, setup_timeout = 30)),
                 error = function(e) NULL)
  if (is.null(cl)) return(NULL)
  ok <- tryCatch({
    parallel::clusterCall(cl, function(lp) {
      .libPaths(lp)
      loadNamespace("psyWiz")
      TRUE
    }, .libPaths())
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    try(parallel::stopCluster(cl), silent = TRUE)
    return(NULL)
  }
  cl
}

#' Run a list of fitting tasks (internal)
#'
#' @param tasks List of task lists, each with a `kind` element.
#' @param parallel `"auto"`, `TRUE` or `FALSE`.
#' @param ncores Maximum number of workers; `NULL` = number of cores minus one.
#' @param threshold Seconds a fit must take before workers are started.
#' @return list(results, mode, workers)
#' @keywords internal
#' @noRd
psw_run_tasks <- function(tasks, parallel = "auto", ncores = NULL, threshold = 1) {
  n <- length(tasks)
  results <- vector("list", n)
  if (identical(parallel, "auto")) mode <- "auto"
  else mode <- if (isTRUE(as.logical(parallel))) "force" else "off"
  done <- 0L

  if (mode == "auto") {
    # fit here until one model turns out to be slow while at least two remain
    slow <- FALSE
    while (done < n && !slow) {
      done <- done + 1L
      results[[done]] <- psw_run_task(tasks[[done]])
      slow <- isTRUE(results[[done]]$time > threshold) && n - done >= 2L
    }
    if (!slow) mode <- "off"
  }

  rest <- setdiff(seq_len(n), seq_len(done))
  workers <- 0L
  if (mode != "off" && length(rest) >= 2L) {
    k <- if (is.null(ncores)) psw_default_cores(length(rest)) else min(as.integer(ncores), length(rest))
    cl <- if (k >= 2L) psw_cluster(k) else NULL
    if (!is.null(cl)) {
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      worker_fun <- function(task) get("psw_run_task", envir = asNamespace("psyWiz"))(task)
      environment(worker_fun) <- globalenv()
      out <- tryCatch(parallel::parLapply(cl, tasks[rest], worker_fun), error = function(e) NULL)
      if (!is.null(out)) {
        results[rest] <- out
        rest <- integer()
        workers <- k
      }
    }
  }

  for (i in rest) results[[i]] <- psw_run_task(tasks[[i]])
  list(results = results, mode = if (workers > 0L) "parallel" else "sequential", workers = workers)
}

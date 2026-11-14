# engine-data.R -- data preparation shared by every psyWiz analysis.
#
# The engine never hands user-supplied names to lavaan. Items become i1..ik and
# factors f1..fm, so names with spaces, accents or operators cannot break the
# model syntax; results are translated back before they are returned.
# Nothing in the engine-*.R files depends on jmvcore.

psw_group_col <- "psw_group"

# Normalise the `factors` argument to a named list of unique item names.
psw_normalise_factors <- function(factors) {
  if (is.character(factors)) factors <- list(Factor = factors)
  if (!is.list(factors) || length(factors) == 0L)
    stop("`factors` must be a named list of item names, e.g. list(F1 = c('x1', 'x2', 'x3')).", call. = FALSE)
  factors <- lapply(factors, function(v) unique(as.character(unlist(v))))
  nms <- names(factors)
  if (is.null(nms)) nms <- rep("", length(factors))
  blank <- is.na(nms) | !nzchar(trimws(nms))
  nms[blank] <- paste("Factor", which(blank))
  names(factors) <- make.unique(trimws(nms), sep = " ")
  factors[lengths(factors) > 0L]
}

# Is a numeric vector made of a handful of whole numbers (i.e. probably Likert)?
psw_looks_ordinal <- function(x, max_levels = 7L) {
  x <- x[!is.na(x)]
  length(x) > 0L && all(x == round(x)) && length(unique(x)) <= max_levels
}

psw_as_numeric <- function(x) {
  if (!is.factor(x)) return(as.numeric(x))
  lv <- suppressWarnings(as.numeric(levels(x)))
  if (anyNA(lv)) lv <- seq_along(levels(x))
  lv[as.integer(x)]
}

psw_as_ordered <- function(x) {
  if (is.factor(x)) return(factor(x, levels = levels(droplevels(x)), ordered = TRUE))
  factor(x, levels = sort(unique(x[!is.na(x)])), ordered = TRUE)
}

#' Prepare data for a psyWiz model (internal)
#'
#' Validates the factor specification, resolves how items are treated
#' (ordinal or continuous), drops unusable rows and renames everything to
#' syntax-safe identifiers.
#'
#' @param data A data frame.
#' @param factors Named list of character vectors of item names.
#' @param group Name of the grouping variable, or `NULL` for a single group.
#' @param type `"auto"`, `"ordinal"` or `"continuous"`.
#' @return A list (`data`, `items`, `fnames`, `type`, `nthr`, `groups`, `n`,
#'   `has_missing`, `notes`).
#' @keywords internal
#' @noRd
psw_prepare <- function(data, factors, group = NULL, type = c("auto", "ordinal", "continuous")) {
  type <- match.arg(type)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  factors <- psw_normalise_factors(factors)
  if (length(factors) == 0L) stop("No items were assigned to any factor.", call. = FALSE)
  notes <- character()

  items <- unlist(factors, use.names = FALSE)
  miss <- setdiff(c(items, group), names(data))
  if (length(miss))
    stop("Variable(s) not found in the data: ", paste(miss, collapse = ", "), call. = FALSE)
  dup <- unique(items[duplicated(items)])
  if (length(dup))
    stop("Each item may load on one factor only (cross-loadings are not supported yet). ",
         "Repeated: ", paste(dup, collapse = ", "), call. = FALSE)
  if (!is.null(group) && group %in% items)
    stop("The grouping variable cannot also be an item.", call. = FALSE)
  if (length(factors) == 1L && length(items) < 3L)
    stop("A one-factor model needs at least 3 items.", call. = FALSE)
  small <- names(factors)[lengths(factors) < 2L]
  if (length(small))
    stop("Every factor needs at least 2 items (3 or more recommended). Too few in: ",
         paste(small, collapse = ", "), call. = FALSE)
  two <- names(factors)[lengths(factors) == 2L]
  if (length(two))
    notes <- c(notes, paste0("Factor(s) with only two items (", paste(two, collapse = ", "),
                             ") are identified only through their correlations with other factors."))

  # ---- how are the items treated? ----
  is_fac <- vapply(data[items], is.factor, logical(1))
  is_num <- vapply(data[items], is.numeric, logical(1))
  bad <- items[!(is_fac | is_num)]
  if (length(bad))
    stop("Items must be numeric or ordinal. Not usable: ", paste(bad, collapse = ", "), call. = FALSE)
  if (type == "auto") {
    if (all(is_fac)) type <- "ordinal"
    else if (all(is_num)) type <- "continuous"
    else stop("The items mix ordinal and continuous variables. Choose whether to treat ",
              "all items as ordinal or all as continuous.", call. = FALSE)
    if (type == "continuous" && all(vapply(data[items], psw_looks_ordinal, logical(1))))
      notes <- c(notes, paste0("All items are stored as continuous but take 7 or fewer whole-number ",
                               "values. If they are Likert-type items, consider treating them as ordinal."))
  }

  # ---- rows ----
  d <- data[items]
  n0 <- nrow(d)
  g <- NULL
  if (!is.null(group)) {
    g <- data[[group]]
    keep <- !is.na(g)
    if (any(!keep))
      notes <- c(notes, sprintf("%d row(s) with a missing value on the grouping variable were dropped.", sum(!keep)))
    d <- d[keep, , drop = FALSE]
    g <- droplevels(as.factor(g[keep]))
  }
  empty <- rowSums(!is.na(d)) == 0L
  if (any(empty)) {
    notes <- c(notes, sprintf("%d row(s) with no item responses were dropped.", sum(empty)))
    d <- d[!empty, , drop = FALSE]
    if (!is.null(g)) g <- droplevels(g[!empty])
  }
  if (nrow(d) == 0L) stop("No usable rows in the data.", call. = FALSE)
  if (!is.null(g) && nlevels(g) < 2L)
    stop("The grouping variable needs at least two groups with data.", call. = FALSE)

  # ---- coerce ----
  nthr <- rep(NA_integer_, length(items))
  if (type == "ordinal") {
    d[] <- lapply(d, psw_as_ordered)
    nlev <- vapply(d, nlevels, integer(1))
    if (any(nlev < 2L))
      stop("Item(s) with a single observed category: ", paste(items[nlev < 2L], collapse = ", "), call. = FALSE)
    if (any(nlev > 10L))
      stop("Item(s) with more than 10 categories cannot be treated as ordinal: ",
           paste(items[nlev > 10L], collapse = ", "), ". Treat the items as continuous.", call. = FALSE)
    nthr <- nlev - 1L
    if (!is.null(g)) {
      # lavaan needs every category of every item to be observed in every group
      for (j in seq_along(items)) {
        tab <- table(g, d[[j]])
        if (any(tab == 0L)) {
          w <- which(tab == 0L, arr.ind = TRUE)[1L, ]
          stop(sprintf(paste0("Item '%s': category '%s' was not observed in group '%s'. Ordinal ",
                              "invariance testing needs every category in every group; merge sparse ",
                              "categories or treat the items as continuous."),
                       items[j], colnames(tab)[w[2L]], rownames(tab)[w[1L]]), call. = FALSE)
        }
      }
    }
  } else {
    d[] <- lapply(d, psw_as_numeric)
    sds <- vapply(d, stats::sd, numeric(1), na.rm = TRUE)
    if (any(is.na(sds) | sds == 0))
      stop("Item(s) with no variance: ", paste(items[is.na(sds) | sds == 0], collapse = ", "), call. = FALSE)
  }

  ids <- paste0("i", seq_along(items))
  fid <- rep(paste0("f", seq_along(factors)), lengths(factors))
  item_map <- data.frame(id = ids, name = items, fid = fid,
                         factor = rep(names(factors), lengths(factors)),
                         stringsAsFactors = FALSE)
  names(d) <- ids
  names(nthr) <- ids
  groups <- NULL
  n <- c(Total = nrow(d))
  if (!is.null(g)) {
    groups <- levels(g)
    n <- table(g)
    n <- stats::setNames(as.integer(n), names(n))
    d[[psw_group_col]] <- as.character(g)
    if (any(n < 100L))
      notes <- c(notes, paste0("Small group(s): ",
                               paste(sprintf("%s (n = %d)", names(n)[n < 100L], n[n < 100L]), collapse = ", "),
                               ". Invariance tests have little power and estimates may be unstable."))
  }
  rownames(d) <- NULL

  list(data = d, items = item_map,
       fnames = stats::setNames(names(factors), paste0("f", seq_along(factors))),
       type = type, nthr = nthr, groups = groups, n = n,
       has_missing = anyNA(d[ids]), n_dropped = n0 - nrow(d), notes = notes)
}

# Tests for the single-scale, ESEM, multiverse and composite engines.
ns <- asNamespace("psyWiz")
data(wizdata, package = "psyWiz")
num <- wizdata
num[-1] <- lapply(num[-1], as.integer)
two <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))

test_that("Guttman and Mokken statistics behave at the known extremes", {
  set.seed(1)
  th <- rnorm(400)
  G <- sapply(c(-1, -0.3, 0.4, 1.1), function(b) as.integer(th > b))
  g <- ns$psw_guttman(G)
  expect_equal(c(g$cr, g$cs, g$errors, g$perfect), c(1, 1, 0, 1))
  expect_equal(ns$psw_mokken_h(G)$H, 1)
  R <- matrix(rbinom(4000, 1, 0.5), 1000)
  expect_lt(abs(ns$psw_mokken_h(R)$H), 0.1)
  expect_lt(ns$psw_guttman(R)$cs, 0.4)
  # binary H is 1 - observed / expected Guttman errors
  x <- G[, 1]; y <- as.integer(th + rnorm(400) > 0.5)
  easy <- if (mean(x) > mean(y)) x else y; hard <- if (mean(x) > mean(y)) y else x
  h <- 1 - mean(hard == 1 & easy == 0) / (mean(hard) * (1 - mean(easy)))
  expect_equal(ns$psw_mokken_h(cbind(x, y))$H, h, tolerance = 1e-10)
})

test_that("reliability coefficients match their definitions", {
  s <- scaleWiz(num, paste0("stress", 1:4), "Stress", estimator = "ML", missing = "listwise", parallel = FALSE)
  X <- as.matrix(stats::na.omit(num[paste0("stress", 1:4)]))
  C <- cov(X)
  a <- 4 / 3 * (1 - sum(diag(C)) / sum(C))
  rel <- setNames(s$reliability$value, s$reliability$coefficient)
  expect_equal(rel[["Cronbach's alpha"]], a)
  expect_gte(rel[["Guttman's lambda-2"]], a)
  fit <- lavaan::cfa("f =~ stress1 + stress2 + stress3 + stress4", num, std.lv = TRUE)
  est <- lavaan::lavInspect(fit, "est")
  expect_equal(rel[["McDonald's omega"]], sum(est$lambda)^2 / (sum(est$lambda)^2 + sum(diag(est$theta))),
               tolerance = 1e-5)
  expect_equal(s$items$alpha.drop[1], 3 / 2 * (1 - sum(diag(C[-1, -1])) / sum(C[-1, -1])))
})

test_that("tau-equivalent and parallel models equal constrained lavaan models", {
  s <- scaleWiz(num, paste0("stress", 1:4), estimator = "ML", missing = "listwise", parallel = FALSE)
  expect_equal(s$fit$step, c("congeneric", "tau", "parallel"))
  tau <- lavaan::cfa("f =~ a*stress1 + a*stress2 + a*stress3 + a*stress4", num, std.lv = TRUE)
  par <- lavaan::cfa(paste("f =~ a*stress1 + a*stress2 + a*stress3 + a*stress4",
                           "stress1 ~~ e*stress1\nstress2 ~~ e*stress2\nstress3 ~~ e*stress3\nstress4 ~~ e*stress4",
                           sep = "\n"), num, std.lv = TRUE)
  expect_equal(s$fit$chisq[2], lavaan::fitMeasures(tau, "chisq")[[1]], tolerance = 1e-4)
  expect_equal(s$fit$chisq[3], lavaan::fitMeasures(par, "chisq")[[1]], tolerance = 1e-4)
  expect_equal(s$fit$df, c(2, 5, 8))
  ord <- scaleWiz(wizdata, paste0("stress", 1:4), parallel = FALSE)
  expect_equal(ord$fit$step, c("congeneric", "tau"))
  expect_true(all(c("Ordinal alpha", "Ordinal omega") %in% ord$reliability$coefficient))
  expect_output(print(ord), "Scalability")
  expect_s3_class(plot(ord), "ggplot")
})

test_that("ESEM reproduces lavaan and compares scoring methods", {
  e <- esemWiz(num, two, estimator = "ML", missing = "listwise", parallel = FALSE)
  expect_equal(e$fit$df, c(19, 13))
  items <- unlist(two)
  ref <- lavaan::cfa(paste0('efa("b")*f1 + efa("b")*f2 =~ ', paste(items, collapse = " + ")), num, rotation = "geomin")
  expect_equal(e$fit$chisq[2], lavaan::fitMeasures(ref, "chisq")[[1]], tolerance = 1e-3)
  own <- e$loadings[e$loadings$target, ]
  expect_true(all(own$esem > 0.5))
  expect_true(all(abs(e$loadings$esem[!e$loadings$target]) < 0.2))
  cc <- e$correlations
  expect_equal(cc$composite, cor(rowMeans(num[two$Trust], na.rm = TRUE), rowMeans(num[two$Stress], na.rm = TRUE)))
  expect_lt(abs(cc$composite), abs(cc$cfa))          # attenuation
  expect_lt(abs(cc$bartlett), abs(cc$scores))        # regression scores over-correlate
  expect_error(esemWiz(num, list(A = two$Trust)), "at least two factors")
  expect_output(print(e), "Factor correlations by method")
})

test_that("the verdict multiverse covers treatments, estimators and rules", {
  mv <- invarianceMultiverse(wizdata, list(Stress = paste0("stress", 1:4)), "country",
                             treatments = "continuous", parallel = FALSE)
  expect_equal(mv$info$specifications, c("Continuous, MLR", "Continuous, ML"))
  expect_equal(nrow(mv$grid), 2 * 5 * 2)
  expect_equal(mv$summary$reading[mv$summary$step == "scalar"], "robustly not supported")
  expect_equal(names(ns$psw_multiverse_wide(mv)), c("specification", "rule", "metric", "scalar"))
})

test_that("factorWiz gives the same CFA as lavaan and the same correlations as esemWiz", {
  f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
  x <- factorWiz(wizdata, f, keep.fits = TRUE)
  expect_s3_class(x, "psw_cfa")
  expect_equal(x$info$estimator, "WLSMV")
  d <- wizdata
  d[unlist(f)] <- lapply(d[unlist(f)], ordered)
  ref <- lavaan::cfa("T =~ trust1 + trust2 + trust3 + trust4\nS =~ stress1 + stress2 + stress3 + stress4",
                     d, std.lv = TRUE, estimator = "WLSMV", missing = "pairwise")
  expect_equal(x$fit$chisq, lavaan::fitMeasures(ref, "chisq.scaled")[[1]], tolerance = 1e-6)
  expect_equal(x$loadings$std, lavaan::standardizedSolution(ref)$est.std[1:8], tolerance = 1e-6)
  e <- esemWiz(wizdata, f, parallel = FALSE)
  expect_equal(x$correlations[c("composite", "corrected", "cfa", "scores")],
               e$correlations[c("composite", "corrected", "cfa", "scores")], tolerance = 1e-6)
  expect_true(all(c("Residual correlation", "Cross-loading") %in% x$suggestions$kind))
  expect_true(all(diff(x$suggestions$mi) <= 0))

  # continuous: omega is the sum-score omega of scaleWiz
  one <- factorWiz(wizdata, f["Trust"], type = "continuous")
  sc <- scaleWiz(wizdata, f$Trust, type = "continuous", parallel = FALSE)
  expect_equal(one$reliability$omega, sc$reliability$value[sc$reliability$coefficient == "McDonald's omega"], tolerance = 1e-6)
  expect_equal(one$reliability$alpha, sc$reliability$value[sc$reliability$coefficient == "Cronbach's alpha"], tolerance = 1e-6)
  expect_null(one$correlations)
  expect_s3_class(plot(one), "ggplot")
  expect_output(print(x), "Factor correlations by scoring method")
})

test_that("the Guttman statistics follow the order of the items, H does not", {
  set.seed(2)
  th <- rnorm(500)
  G <- as.data.frame(sapply(c(-1.2, -0.4, 0.3, 1.1), function(b) as.integer(th + rnorm(500, 0, 0.3) > b)))
  names(G) <- paste0("g", 1:4)
  right <- scaleWiz(G, paste0("g", 1:4), type = "continuous", parallel = FALSE)
  wrong <- scaleWiz(G, paste0("g", 4:1), type = "continuous", parallel = FALSE)
  cr <- function(x, order) x$guttman[[order]]$table$value[4]
  # entered easiest first = the sorted order; entered backwards is much worse
  expect_equal(cr(right, "entered"), cr(right, "difficulty"))
  expect_lt(cr(wrong, "entered"), cr(right, "entered") - 0.2)
  expect_equal(cr(wrong, "difficulty"), cr(right, "difficulty"))
  tau <- function(x) utils::tail(x$scalability$value, 1)
  expect_equal(c(tau(right), tau(wrong)), c(1, -1))
  expect_equal(right$scalability$value[1], wrong$scalability$value[1])
  expect_s3_class(plot(right, type = "scalogram"), "ggplot")
  expect_true(right$guttman$entered$scalogram$binary)
})

test_that("CFA vs ESEM reports the uncorrected chi-square, which cannot favour the CFA", {
  e <- esemWiz(wizdata, two, parallel = FALSE)
  expect_lte(e$fit$chisq.std[2], e$fit$chisq.std[1])
  expect_false(e$info$reversed)
})

test_that("plots omit internal instructional captions and attach model fit cleanly", {
  f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
  cfa <- factorWiz(wizdata, f)
  p_plain <- plot(cfa, type = "comb")
  expect_null(p_plain$labels$caption)
  p_below <- plot(cfa, type = "comb", fit = "below")
  expect_match(p_below$labels$caption, "Estimator: WLSMV", fixed = TRUE)
  p_above <- plot(cfa, type = "comb", fit = "above")
  expect_match(p_above$labels$subtitle, "Estimator: WLSMV", fixed = TRUE)
})

test_that("bifactorWiz recovers a known structure and its indices follow from the loadings", {
  set.seed(5)
  n <- 3000
  g <- rnorm(n)
  d <- do.call(cbind, lapply(1:3, function(k) { s <- rnorm(n); sapply(1:4, function(j) 0.6 * g + 0.45 * s + rnorm(n, 0, 0.66)) }))
  colnames(d) <- paste0(rep(c("a", "b", "c"), each = 4), 1:4)
  f <- list(A = paste0("a", 1:4), B = paste0("b", 1:4), C = paste0("c", 1:4))
  x <- bifactorWiz(as.data.frame(d), f)
  # population: variance = .36 + .2025 + .4356 = .998, so loadings are about .60 and .45
  expect_equal(mean(x$loadings$general), 0.60, tolerance = 0.03)
  expect_equal(mean(x$loadings$own), 0.45, tolerance = 0.04)
  L <- x$loadings
  ss <- tapply(L$own, L$specific, sum)
  V <- sum(L$general)^2 + sum(ss^2) + sum(L$unique)
  expect_equal(x$indices$omegaH[1], sum(L$general)^2 / V)
  expect_equal(x$indices$ecv[1], sum(L$general^2) / sum(L$general^2 + L$own^2))
  expect_equal(x$info$puc, 1 - 3 * 6 / 66)
  expect_lt(x$fit$chisq[1], x$fit$chisq[2])
  for (t in c("antenna", "facing", "path")) expect_s3_class(plot(x, type = t), "ggplot")
  expect_output(print(x), "omega hierarchical")
  expect_error(bifactorWiz(wizdata, list(Trust = paste0("trust", 1:2), Stress = paste0("stress", 1:4))), "at least three items")

  # second-order model on the same data: equivalent to three correlated factors, so same fit
  h <- factorWiz(as.data.frame(d), f, higher = "G", residuals = list(c("a1", "a2")))
  cf <- factorWiz(as.data.frame(d), f, residuals = list(c("a1", "a2")))
  expect_equal(h$fit$chisq, cf$fit$chisq, tolerance = 1e-4)
  expect_equal(h$correlations$cfa, cf$correlations$cfa, tolerance = 1e-3)
  expect_equal(h$reliability$omega, cf$reliability$omega, tolerance = 1e-3)
  expect_equal(nrow(h$residuals), 1)
  expect_s3_class(plot(h, thick = TRUE), "ggplot")
  expect_s3_class(plot(cf, type = "path", numbers = FALSE), "ggplot")
})

test_that("the brush plot is drawn from an ESEM", {
  e <- esemWiz(wizdata, two, parallel = FALSE)
  expect_equal(names(e$uniqueness), unlist(two, use.names = FALSE))
  expect_true(all(e$uniqueness > 0 & e$uniqueness < 1))
  expect_s3_class(plot(e), "ggplot")
})

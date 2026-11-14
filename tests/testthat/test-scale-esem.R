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

test_that("the composite check accepts one channel and rejects two", {
  set.seed(3)
  n <- 600
  X <- matrix(rnorm(n * 4), n, dimnames = list(NULL, paste0("x", 1:4)))
  comp <- X %*% c(0.6, 0.5, 0.3, 0.2)
  one <- data.frame(X, y1 = 0.5 * comp + rnorm(n), y2 = 0.3 * comp + rnorm(n), y3 = -0.4 * comp + rnorm(n))
  two_ch <- data.frame(X, y1 = 0.5 * X[, 1] + 0.4 * X[, 2] + rnorm(n), y2 = 0.5 * X[, 3] + 0.4 * X[, 4] + rnorm(n))
  a <- compositeWiz(one, paste0("x", 1:4), paste0("y", 1:3))
  b <- compositeWiz(two_ch, paste0("x", 1:4), c("y1", "y2"))
  expect_gt(a$test$p, 0.05)
  expect_lt(b$test$p, 0.001)
  expect_equal(a$test$df, 6)
  expect_true(all(a$r2$retained > 0.9))
  expect_true(all(a$r2$free >= a$r2$composite - 1e-12))
  # in-sample R2 of the free model is the usual regression R2
  expect_equal(a$r2$free[1], summary(lm(y1 ~ x1 + x2 + x3 + x4, one))$r.squared)
  expect_identical(compositeWiz(one, paste0("x", 1:4), paste0("y", 1:3))$r2, a$r2)   # seeded
  expect_error(compositeWiz(one, paste0("x", 1:4), "y1"), "At least two outcomes")
})

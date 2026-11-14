# Engine tests: results are checked against plain lavaan and against known
# properties of the models. Run with testthat::test_local() or R CMD check.
ns <- asNamespace("psyWiz")
data(wizdata, package = "psyWiz")
two <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))
num <- wizdata
num[-1] <- lapply(num[-1], as.integer)

test_that("data preparation validates and resolves the item type", {
  p <- ns$psw_prepare(wizdata, two, "country")
  expect_equal(p$type, "ordinal")
  expect_equal(unname(p$nthr), rep(4L, 8))
  expect_equal(p$groups, c("Chile", "Brazil", "Spain"))
  expect_equal(ns$psw_prepare(num, two, "country")$type, "continuous")
  mixed <- wizdata; mixed$trust1 <- as.integer(mixed$trust1)
  expect_error(ns$psw_prepare(mixed, two, "country"), "mix ordinal and continuous")
  expect_error(ns$psw_prepare(wizdata, list(A = c("trust1", "trust2")), "country"), "at least 3 items")
  expect_error(ns$psw_prepare(wizdata, list(A = paste0("trust", 1:3), B = c("trust1", "stress1")), "country"),
               "one factor only")
  expect_error(ns$psw_prepare(wizdata, list(A = c("trust1", "nope", "trust2"))), "not found")
  sparse <- wizdata[!(wizdata$country == "Spain" & wizdata$trust1 == 5 & !is.na(wizdata$trust1)), ]
  expect_error(ns$psw_prepare(sparse, two, "country"), "not observed in group 'Spain'")
})

test_that("continuous ladder reproduces lavaan's group.equal models", {
  res <- invarianceWiz(num, two, "country", estimator = "ML", missing = "listwise",
                       strict = TRUE, parallel = FALSE)
  model <- "Trust =~ trust1 + trust2 + trust3 + trust4\nStress =~ stress1 + stress2 + stress3 + stress4"
  eq <- list(character(), "loadings", c("loadings", "intercepts"), c("loadings", "intercepts", "residuals"))
  for (i in seq_along(eq)) {
    ref <- lavaan::cfa(model, num, group = "country", group.equal = eq[[i]])
    fm <- lavaan::fitMeasures(ref, c("chisq", "df", "cfi", "rmsea"))
    expect_equal(res$fit$chisq[i], fm[["chisq"]], tolerance = 1e-4)
    expect_equal(res$fit$df[i], fm[["df"]])
    expect_equal(res$fit$cfi[i], fm[["cfi"]], tolerance = 1e-5)
  }
  expect_equal(res$fit$step, c("configural", "metric", "scalar", "strict"))
  expect_equal(res$comparison$verdict, c("Supported", "Not supported", "Supported"))
  # the simulated violation is the intercept of stress3 in Brazil
  top <- res$culprits[res$culprits$step == "scalar", ][1, ]
  expect_equal(c(top$item, top$group), c("stress3", "Brazil"))
})

test_that("ordinal ladder has the degrees of freedom of Wu & Estabrook (2016)", {
  res <- invarianceWiz(wizdata, two, "country", parallel = FALSE)
  expect_equal(res$info$estimator, "WLSMV")
  expect_equal(res$fit$step, c("configural", "thresholds", "metric", "scalar"))
  # p = 8 items, T = 4 thresholds, m = 2 factors, G = 3 groups
  expect_equal(res$fit$df, c(57, 57 + 8 * 2 * 2, 89 + 6 * 2, 101 + 6 * 2))
  expect_equal(res$comparison$verdict, c("Supported", "Supported", "Not supported"))
  top <- res$culprits[res$culprits$step == "scalar", ][1, ]
  expect_equal(c(top$item, top$group), c("stress3", "Brazil"))
  expect_true(all(res$info$variant == "scaled"))
})

test_that("with three categories the thresholds model equals the configural model", {
  d3 <- wizdata
  d3[-1] <- lapply(d3[-1], function(x) factor(pmin(pmax(as.integer(x), 2L), 4L), ordered = TRUE))
  res <- invarianceWiz(d3, list(Stress = paste0("stress", 1:4)), "country", parallel = FALSE)
  expect_equal(res$fit$df[1], res$fit$df[2])
  expect_equal(res$fit$chisq[1], res$fit$chisq[2], tolerance = 1e-3)
  expect_equal(res$comparison$verdict[1], "Equivalent to previous model")
})

test_that("verdicts follow the criteria and cutoffs without refitting", {
  res <- invarianceWiz(num, two, "country", parallel = FALSE)
  strictest <- invarianceVerdict(res, criteria = "cfi", cutoffs = list(cfi = 0))
  expect_true(all(strictest$comparison$flag.cfi == (strictest$comparison$dcfi < 0)))
  lenient <- invarianceVerdict(res, criteria = c("cfi", "rmsea"), cutoffs = list(cfi = 0.5, rmsea = 0.5))
  expect_true(all(lenient$comparison$verdict == "Supported"))
  mixed <- invarianceVerdict(res, criteria = c("cfi", "rmsea"), cutoffs = list(cfi = 0.5, rmsea = 0.015))
  expect_equal(mixed$comparison$verdict[2], "Ambiguous")
  expect_equal(invarianceVerdict(res, criteria = character())$comparison$verdict[1], "No criterion selected")
})

test_that("multi-estimator results equal separate lavaan fits", {
  res <- cfaWiz(wizdata, two, estimators = c("ML", "MLR", "WLSMV", "ULS"), missing = "listwise", parallel = FALSE)
  model <- "Trust =~ trust1 + trust2 + trust3 + trust4\nStress =~ stress1 + stress2 + stress3 + stress4"
  expect_equal(res$fit$treated, c("continuous", "continuous", "ordinal", "ordinal"))
  ml <- lavaan::cfa(model, num, std.lv = TRUE)
  expect_equal(res$fit$chisq[1], lavaan::fitMeasures(ml, "chisq")[[1]], tolerance = 1e-5)
  wl <- lavaan::cfa(model, wizdata, std.lv = TRUE, estimator = "WLSMV", ordered = unlist(two))
  expect_equal(res$fit$chisq[3], lavaan::fitMeasures(wl, "chisq.scaled")[[1]], tolerance = 1e-5)
  expect_equal(res$fit$cfi[3], lavaan::fitMeasures(wl, "cfi.scaled")[[1]], tolerance = 1e-6)
  std <- lavaan::standardizedSolution(wl)
  expect_equal(res$loadings$std[res$loadings$estimator == "WLSMV"], std$est.std[std$op == "=~"], tolerance = 1e-6)
  rel <- res$reliability[res$reliability$estimator == "WLSMV" & res$reliability$factor == "Trust", ]
  l <- std$est.std[std$op == "=~" & std$lhs == "Trust"]
  expect_equal(rel$omega, sum(l)^2 / (sum(l)^2 + sum(1 - l^2)))
  expect_equal(rel$ave, mean(l^2))
  expect_error(cfaWiz(wizdata, two, estimators = "BOGUS"), "Unknown estimator")
})

test_that("awkward variable and factor names survive", {
  d <- wizdata
  names(d) <- c("país (grupo)", paste("ítem", 1:4, "- confianza"), paste0("i", 1:4))
  f <- list("Confianza social" = names(d)[2:5], "f1" = names(d)[6:9])
  res <- invarianceWiz(d, f, "país (grupo)", parallel = FALSE)
  expect_true(all(res$fit$converged))
  expect_setequal(unique(res$loadings$item), names(d)[2:9])
  expect_setequal(unique(res$loadings$factor), names(f))
  expect_match(res$syntax$configural, "Confianza social =~", fixed = TRUE)
  expect_output(print(res, colored = FALSE), "Intercepts \\(scalar\\)")
})

test_that("parallel and sequential fitting agree", {
  skip_on_cran()
  a <- invarianceWiz(wizdata, two, "country", parallel = FALSE)
  b <- invarianceWiz(wizdata, two, "country", parallel = TRUE, ncores = 2)
  skip_if(b$info$mode != "parallel", "worker processes could not be started")
  expect_equal(a$fit$chisq, b$fit$chisq)
  expect_equal(a$comparison$dchisq, b$comparison$dchisq)
  expect_equal(a$culprits$x2, b$culprits$x2)
})

test_that("fast score tests equal lavaan::lavTestScore", {
  prep <- ns$psw_prepare(wizdata, two, "country")
  fit <- ns$psw_fit(ns$psw_syntax(prep, "metric"), prep, "WLSMV", "pairwise")$fit
  ref <- suppressWarnings(lavaan::lavTestScore(fit))$uni$X2
  expect_equal(ns$psw_score_uni(fit, seq_along(ref)), as.numeric(ref), tolerance = 1e-6)
})

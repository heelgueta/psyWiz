# The jamovi analyses, run headless through the functions jmvtools generates.
data(wizdata, package = "psyWiz")
fl <- list(list(label = "Trust", vars = paste0("trust", 1:4)),
           list(label = "Stress", vars = paste0("stress", 1:4)))

test_that("nothing is fitted until Run is ticked", {
  r <- groupCfaCorr(wizdata, factors = fl, group = "country")
  expect_equal(r$fit$rowCount, 0)
  expect_equal(r$groups$asDF$n, c(420, 380, 400))
  expect_match(r$todo$content, "Tick <b>Run analysis</b>", fixed = TRUE)
  expect_equal(multiEstCfa(wizdata, factors = fl)$fit$rowCount, 0)
})

test_that("group analyses fill their tables", {
  r <- groupCfaCorr(wizdata, factors = fl, group = "country", run = TRUE, parallel = FALSE,
                    showLoadings = TRUE, showDiagram = TRUE, highlight = TRUE)
  expect_equal(r$fit$asDF$model, c("Configural", "Thresholds", "Loadings (metric)", "Intercepts (scalar)"))
  expect_equal(r$comparison$asDF$verdict, c("Supported", "Supported", "Not supported"))
  expect_equal(r$culprits$asDF$item[1], "stress3")
  expect_equal(dim(r$loadings$asDF), c(8L, 5L))
  expect_equal(nrow(r$diagram$state$loadings), 24)

  one <- groupCfaOne(wizdata, vars = paste0("trust", 1:4), group = "country", run = TRUE,
                     parallel = FALSE, type = "continuous", strict = TRUE, culprits = "all", culpritsN = 1)
  expect_equal(nrow(one$fit$asDF), 4)
  expect_equal(nrow(one$culprits$asDF), 3)
})

test_that("multi-estimator analysis fills one column per estimator", {
  r <- multiEstCfa(wizdata, factors = fl, run = TRUE, parallel = FALSE, estML = TRUE, estULSMV = FALSE)
  expect_equal(r$fit$asDF$estimator, c("ML", "MLR", "WLSMV"))
  expect_equal(names(r$std$asDF), c("factor", "item", "ML", "MLR", "WLSMV"))
  expect_equal(nrow(r$cors$asDF), 1)
  expect_equal(ncol(r$rel$asDF), 7)
})

test_that("engine errors are shown as jamovi errors", {
  bad <- list(list(label = "A", vars = paste0("trust", 1:3)), list(label = "B", vars = c("trust1", "stress1")))
  expect_error(groupCfaCorr(wizdata, factors = bad, group = "country", run = TRUE), "one factor only")
})

test_that("scale check, CFA vs ESEM and composite check fill their tables", {
  expect_equal(scaleCheck(wizdata, vars = paste0("trust", 1:4))$fit$rowCount, 0)
  s <- scaleCheck(wizdata, vars = paste0("trust", 1:4), scaleName = "Trust", run = TRUE, parallel = FALSE,
                  showDiagram = TRUE, highlight = TRUE)
  expect_equal(s$fit$asDF$model, c("Congeneric", "Tau-equivalent"))
  expect_equal(nrow(s$items$asDF), 4)
  expect_equal(nrow(s$scalability$asDF), 7)
  expect_equal(s$diagram$state$type, "comb")

  e <- esemCfa(wizdata, factors = fl, run = TRUE, parallel = FALSE)
  expect_equal(nrow(e$fit$asDF), 2)
  expect_equal(names(e$loadings$asDF), c("own", "item", "cfa", "e1", "e2"))
  expect_equal(nrow(e$cors$asDF), 1)

  d <- data.frame(matrix(rnorm(300 * 5), 300))
  cc <- compositeCheck(d, indicators = c("X1", "X2", "X3"), outcomes = c("X4", "X5"))
  expect_equal(nrow(cc$r2$asDF), 2)
  expect_equal(nrow(cc$weights$asDF), 3)
})

test_that("the robustness tables are filled on request", {
  g <- groupCfaOne(wizdata, vars = paste0("stress", 1:4), group = "country", run = TRUE, parallel = FALSE,
                   type = "continuous", robust = TRUE, showDiagram = TRUE, diagramStyle = "tiles")
  expect_equal(nrow(g$robustSummary$asDF), 3)
  expect_equal(nrow(g$robustGrid$asDF), 20)
  expect_equal(g$diagram$state$type, "tiles")
})

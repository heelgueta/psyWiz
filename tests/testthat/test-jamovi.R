# The jamovi analyses, run headless through the functions jmvtools generates.
data(wizdata, package = "psyWiz")
fl <- list(list(label = "Trust", vars = paste0("trust", 1:4)),
           list(label = "Stress", vars = paste0("stress", 1:4)))

test_that("nothing is fitted until Run is ticked", {
  r <- groupCfa(wizdata, factors = fl, group = "country")
  expect_equal(r$fit$rowCount, 0)
  expect_equal(r$groups$asDF$n, c(420, 380, 400))
  expect_match(r$todo$content, "Tick <b>Run analysis</b>", fixed = TRUE)
  expect_equal(multiEstCfa(wizdata, factors = fl)$fit$rowCount, 0)
})

test_that("group analyses fill their tables", {
  r <- groupCfa(wizdata, factors = fl, group = "country", run = TRUE, parallel = FALSE,
                    showLoadings = TRUE, showComb = TRUE, showPath = TRUE, highlight = TRUE)
  expect_equal(r$fit$asDF$model, c("Configural", "Thresholds", "Loadings (metric)", "Intercepts (scalar)"))
  expect_equal(r$comparison$asDF$verdict, c("Supported", "Supported", "Not supported"))
  expect_equal(r$culprits$asDF$item[1], "stress3")
  expect_equal(dim(r$loadings$asDF), c(8L, 5L))
  expect_equal(nrow(r$pathFig$state$args$loadings), 24)
  expect_equal(r$comb$state$args$type, "comb")

  one <- groupCfa(wizdata, factors = fl[1], group = "country", run = TRUE,
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
  h <- multiEstCfa(wizdata, factors = fl, run = TRUE, parallel = FALSE, highlight = TRUE)
  expect_equal(nrow(h$fit$asDF), 3)
  expect_false("heat" %in% h$itemNames)
})

test_that("engine errors are shown as jamovi errors", {
  bad <- list(list(label = "A", vars = paste0("trust", 1:3)), list(label = "B", vars = c("trust1", "stress1")))
  expect_error(groupCfa(wizdata, factors = bad, group = "country", run = TRUE), "one factor only")
})

test_that("scale check and CFA vs ESEM fill their tables and plot guides", {
  expect_equal(scaleCheck(wizdata, vars = paste0("trust", 1:4))$reliability$rowCount, 0)
  # defaults: reliability, items, and comb plot + guide
  s0 <- scaleCheck(wizdata, vars = paste0("trust", 1:4), run = TRUE, parallel = FALSE)
  expect_gte(nrow(s0$reliability$asDF), 4)
  expect_equal(nrow(s0$itemStats$asDF), 4)
  expect_equal(s0$fit$rowCount, 0)
  expect_equal(s0$scalability$rowCount, 0)
  expect_equal(s0$comb$state$kind, "figure")
  expect_match(s0$combDesc$content, "Reading the Comb Plot", fixed = TRUE)
  expect_lt(match("combDesc", s0$itemNames), match("comb", s0$itemNames))
  s <- scaleCheck(wizdata, vars = paste0("trust", 1:4), scaleName = "Trust", run = TRUE, parallel = FALSE,
                  tau = TRUE, showScal = TRUE, showComb = TRUE, combThick = TRUE, highlight = TRUE, showPlotGuide = FALSE)
  expect_equal(s$fit$asDF$model, c("Congeneric", "Tau-equivalent"))
  expect_equal(nrow(s$itemStats$asDF), 4)
  expect_equal(nrow(s$scalability$asDF), 8)   # 7 statistics + agreement of the entered order
  expect_false(s$combDesc$visible)
  s2 <- scaleCheck(wizdata, vars = paste0("trust", 1:4), run = TRUE, parallel = FALSE, showScal = TRUE,
                   guttmanOrder = "difficulty", showScalogram = TRUE)
  expect_equal(nrow(s2$scalability$asDF), 7)
  expect_equal(dim(s2$scalogram$state$args$sc$z), c(150L, 4L))
  expect_match(s2$scalogramDesc$content, "Reading the Scalogram", fixed = TRUE)
  expect_lt(match("scalogramDesc", s2$itemNames), match("scalogram", s2$itemNames))
  expect_true(s$comb$state$args$opts$thick)

  e <- esemCfa(wizdata, factors = fl, run = TRUE, parallel = FALSE, showBrush = TRUE)
  expect_equal(e$brush$state$kind, "brush")
  expect_match(e$brushDesc$content, "Reading the Brush Plot", fixed = TRUE)
  expect_lt(match("brushDesc", e$itemNames), match("brush", e$itemNames))
  expect_equal(nrow(e$fit$asDF), 2)
  expect_equal(names(e$loadings$asDF), c("own", "item", "cfa", "e1", "e2"))
  expect_equal(nrow(e$cors$asDF), 1)
})

test_that("the robustness tables are filled on request", {
  g <- groupCfa(wizdata, factors = fl[2], group = "country", run = TRUE, parallel = FALSE,
                type = "continuous", robust = TRUE, partial = TRUE, showPath = TRUE, pathNumbers = FALSE)
  expect_equal(g$partialLadder$asDF$verdict, c("Supported", "Supported"))
  # the loadings step is ambiguous for this factor, so the loading of stress3 goes first and its intercept follows
  expect_equal(g$partialReleased$asDF$item, c("stress3", "stress3"))
  expect_null(g$fit$state$prep)
  expect_equal(nrow(g$robustSummary$asDF), 3)
  expect_equal(nrow(g$robustGrid$asDF), 20)
  expect_false(g$pathFig$state$args$opts$numbers)
})

test_that("plain CFA fills its tables, for one factor and for several", {
  expect_equal(plainCfa(wizdata, factors = fl, resCov = list())$fit$rowCount, 0)
  r <- plainCfa(wizdata, factors = fl, resCov = list(), run = TRUE, showMethods = TRUE, showMi = TRUE, miN = 3, highlight = TRUE)
  expect_equal(nrow(r$fit$asDF), 1)
  expect_equal(nrow(r$loadings$asDF), 8)
  expect_equal(nrow(r$cors$asDF), 1)
  expect_equal(r$rel$asDF$k, c(4, 4))
  expect_equal(nrow(r$methods$asDF), 1)
  expect_equal(nrow(r$mi$asDF), 3)
  expect_equal(r$comb$state$kind, "figure")   # the comb plot is on by default
  expect_match(r$combDesc$content, "Reading the Comb Plot", fixed = TRUE)
  expect_lt(match("combDesc", r$itemNames), match("comb", r$itemNames))
  expect_false("heat" %in% r$itemNames)

  one <- plainCfa(wizdata, factors = fl[1], resCov = list(), run = TRUE, type = "continuous", showMethods = TRUE)
  expect_equal(nrow(one$loadings$asDF), 4)
  expect_false(one$cors$visible)
  expect_false(one$methods$visible)
})

# (called from R, the Pairs option `resCov` has no default, as in jmv::cfa: pass list())

# Simulated questionnaire with a general factor and three specific ones, used for the
# second-order and bifactor analyses (the example data have only two factors).
bf_data <- function(n = 600) {
  set.seed(5)
  g <- rnorm(n)
  d <- do.call(cbind, lapply(1:3, function(k) {
    s <- rnorm(n)
    sapply(1:4, function(j) 0.6 * g + 0.45 * s + rnorm(n, 0, 0.66))
  }))
  colnames(d) <- paste0(rep(c("a", "b", "c"), each = 4), 1:4)
  as.data.frame(d)
}
bfl <- lapply(c("a", "b", "c"), function(k) list(label = toupper(k), vars = paste0(k, 1:4)))

test_that("plain CFA takes residual correlations and a second-order factor", {
  d <- bf_data()
  r <- plainCfa(d, factors = bfl, run = TRUE, secondOrder = TRUE, higherName = "G",
                resCov = list(list(i1 = "a1", i2 = "a2")), showMethods = TRUE)
  expect_equal(r$higher$asDF$factor, c("A", "B", "C"))
  expect_true(all(r$higher$asDF$std > 0.5))
  expect_equal(unlist(r$resid$asDF[1, 1:2], use.names = FALSE), c("a1", "a2"))
  expect_false(r$cors$visible)
  expect_equal(nrow(r$methods$asDF), 3)
  expect_equal(r$comb$state$args$opts$higher$name, "G")
  expect_equal(nrow(r$comb$state$args$opts$resid), 1)
  expect_error(plainCfa(wizdata, factors = fl, resCov = list(), run = TRUE, secondOrder = TRUE), "at least three factors")
})

test_that("bifactor CFA fills its tables and figures", {
  d <- bf_data()
  expect_equal(bifactorCfa(d, factors = bfl)$fit$rowCount, 0)
  b <- bifactorCfa(d, factors = bfl, run = TRUE, showFacing = TRUE, showPath = TRUE, highlight = TRUE)
  expect_equal(b$fit$asDF$model, c("Bifactor", "Correlated factors"))
  expect_equal(b$indices$asDF$scale, c("Total score", "A", "B", "C"))
  expect_equal(nrow(b$loadings$asDF), 12)
  expect_equal(c(b$antenna$state$kind, b$facing$state$kind, b$pathFig$state$kind), c("antenna", "facing", "bfpath"))
  expect_match(b$antennaDesc$content, "Reading the Antenna Plot", fixed = TRUE)
  expect_lt(match("antennaDesc", b$itemNames), match("antenna", b$itemNames))
  expect_match(b$facingDesc$content, "Reading Facing Combs", fixed = TRUE)
  expect_lt(match("facingDesc", b$itemNames), match("facing", b$itemNames))
})

test_that("no result is named like a member of jamovi's results classes", {
  # such a name replaces jamovi's own member and breaks the analysis in the jamovi window only
  reserved <- c(names(jmvcore::Group$active), names(jmvcore::Group$public_methods),
                names(jmvcore::Group$get_inherit()$active), names(jmvcore::Group$get_inherit()$public_methods))
  for (r in list(scaleCheck(wizdata, vars = paste0("trust", 1:4)),
                 plainCfa(wizdata, factors = fl, resCov = list()), groupCfa(wizdata, factors = fl, group = "country"),
                 bifactorCfa(wizdata, factors = fl), esemCfa(wizdata, factors = fl), multiEstCfa(wizdata, factors = fl)))
    expect_length(intersect(r$itemNames, reserved), 0)
})

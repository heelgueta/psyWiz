data(wizdata, package = "psyWiz")
fl <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))

test_that("invarianceMeans handles unsupported scalar invariance", {
  res <- invarianceWiz(wizdata, fl, group = "country", parallel = FALSE)
  # Full scalar fails on wizdata without partial releases
  expect_equal(res$comparison$verdict[res$comparison$step == "scalar"], "Not supported")
  m <- invarianceMeans(res)
  expect_false(attr(m, "supported"))
  expect_equal(nrow(m), 0)
  expect_match(attr(m, "notes"), "scalar invariance is not supported", fixed = TRUE)

  # Plotting unsupported returns informative message plot
  p <- plot(m)
  expect_s3_class(p, "ggplot")
  expect_equal(attr(p, "psw_size"), c(width = 520, height = 240))
})

test_that("invarianceMeans computes full scalar means when supported", {
  # wizdata continuous with Trust factor only: scalar is supported
  res <- invarianceWiz(wizdata, list(Trust = paste0("trust", 1:4)), group = "country",
                       type = "continuous", parallel = FALSE, means = TRUE)
  expect_true("scalar" %in% res$comparison$step)
  expect_equal(res$comparison$verdict[res$comparison$step == "scalar"], "Supported")

  m <- res$means
  expect_true(attr(m, "supported"))
  expect_equal(attr(m, "model_type"), "full")
  expect_equal(attr(m, "reference"), "Chile")
  expect_equal(nrow(m), 3) # Chile, Brazil, Spain

  chile <- m[m$group == "Chile", ]
  expect_true(chile$is_ref)
  expect_equal(chile$mean, 0)
  expect_true(is.na(chile$se))
  expect_true(is.na(chile$z))
  expect_true(is.na(chile$p))
  expect_equal(chile$d_ref, 0)
  expect_equal(chile$d_pool, 0)

  brazil <- m[m$group == "Brazil", ]
  expect_false(brazil$is_ref)
  expect_lt(brazil$mean, 0)
  expect_gt(brazil$se, 0)
  expect_lt(brazil$z, 0)
  expect_lt(brazil$p, 0.05)
  expect_equal(brazil$d_ref, brazil$mean, tolerance = 1e-4) # ref SD is 1.0

  # Custom reference group
  m_spain <- invarianceMeans(res, reference = "Spain")
  expect_equal(attr(m_spain, "reference"), "Spain")
  spain <- m_spain[m_spain$group == "Spain", ]
  expect_true(spain$is_ref)
  expect_equal(spain$mean, 0)

  # Invalid reference group
  expect_error(invarianceMeans(res, reference = "Nonexistent"), "not found")

  # Plot produces ggplot with psw_size attribute
  p <- plot(m)
  expect_s3_class(p, "ggplot")
  expect_true(!is.null(attr(p, "psw_size")))
})

test_that("invarianceMeans supports partial scalar invariance", {
  # On 2-factor wizdata, full scalar fails but partial scalar holds when stress3 is released
  res <- invarianceWiz(wizdata, fl, group = "country", parallel = FALSE,
                       partial = TRUE, means = TRUE)
  m <- res$means
  expect_true(attr(m, "supported"))
  expect_equal(attr(m, "model_type"), "partial")
  expect_equal(attr(m, "freed_items"), "stress3")
  expect_equal(nrow(m), 6) # 2 factors * 3 groups

  # Output printed cleanly
  expect_output(print(m), "partial scalar invariance")
  expect_output(print(m), "stress3")

  # Directly passing partial object
  part <- res$partial
  m_part <- invarianceMeans(part)
  expect_true(attr(m_part, "supported"))
  expect_equal(attr(m_part, "model_type"), "partial")
})

test_that("invarianceMeans works with ordinal item treatment", {
  res <- invarianceWiz(wizdata, list(Trust = paste0("trust", 1:4)), group = "country",
                       type = "ordinal", parallel = FALSE, means = TRUE)
  expect_equal(res$comparison$verdict[res$comparison$step == "scalar"], "Supported")
  m <- res$means
  expect_true(attr(m, "supported"))
  expect_equal(nrow(m), 3)
  expect_equal(m$mean[m$group == "Chile"], 0)
})

test_that("groupCfa jamovi analysis populates latent means and plot", {
  fl_jmv <- list(list(label = "Trust", vars = paste0("trust", 1:4)))
  r <- groupCfa(wizdata, factors = fl_jmv, group = "country", run = TRUE, parallel = FALSE,
                type = "continuous", latentMeans = TRUE, latentMeansPlot = TRUE)
  df <- r$latentMeans$asDF
  expect_equal(nrow(df), 3)
  expect_equal(df$group, c("Chile", "Brazil", "Spain"))
  expect_equal(df$mean[1], 0)
  expect_true(is.na(df$se[1]))
  expect_false(is.na(df$se[2]))
  expect_equal(r$latentMeansPlot$state$kind, "means")
})

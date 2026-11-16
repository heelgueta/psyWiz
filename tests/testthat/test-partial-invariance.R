# Partial invariance: the example data have one non-invariant intercept
# (stress3 in Brazil), which the search has to find and release.
data(wizdata, package = "psyWiz")
f <- list(Trust = paste0("trust", 1:4), Stress = paste0("stress", 1:4))

test_that("the search releases the non-invariant intercept, continuous and ordinal", {
  res_c <- invarianceWiz(wizdata, f, group = "country", type = "continuous", partial = TRUE, parallel = FALSE)
  part <- res_c$partial
  expect_s3_class(part, "psw_partial")
  expect_equal(part$freed, list(load = character(), int = "stress3", res = character()))
  expect_equal(part$comparison$verdict, c("Supported", "Supported"))
  expect_equal(part$comparison$freed, c("", "stress3"))
  expect_match(part$summary, "Partial scalar invariance is supported")
  expect_equal(part$released$group, "Brazil")
  # the steps before the release are the same models as in the full ladder
  expect_equal(part$fit$chisq[1:2], res_c$fit$chisq[1:2])

  res_o <- invarianceWiz(wizdata, f, group = "country", type = "ordinal", parallel = FALSE)
  part_o <- invariancePartial(res_o)
  expect_equal(part_o$freed$int, "stress3")
  expect_equal(part_o$comparison$verdict[part_o$comparison$step == "scalar"], "Supported")
  expect_true(any(grepl("partial measurement invariance", capture.output(print(part_o)))))

  # releasing by hand gives the same model as the search
  man <- invariancePartial(res_c, auto = FALSE, free = list(int = "stress3"))
  expect_equal(man$fit$chisq, part$fit$chisq)
  expect_equal(man$released$reason, "Chosen by the user")
  expect_error(invariancePartial(res_c, free = list(int = "nope")), "Not an item")
  expect_error(invariancePartial(res_c, free = list(thr = "stress3")), "may only contain")
})

test_that("limits on releases are respected and a freed loading frees the intercept", {
  # three items per factor and two must stay constrained: at most one release per factor
  f3 <- list(Trust = paste0("trust", 1:3), Stress = paste0("stress", 1:3))
  res <- invarianceWiz(wizdata, f3, group = "country", type = "continuous", parallel = FALSE)
  # an impossible cutoff makes every step fail, so the search goes as far as it is allowed to
  part <- invariancePartial(res, max.free = 10L, cutoffs = list(cfi = -1, rmsea = -1))
  per_factor <- function(items) table(factor(sub("[0-9]$", "", items), levels = c("trust", "stress")))
  expect_true(all(per_factor(part$freed$load) == 1L))
  expect_true(all(part$freed$load %in% part$freed$int))
  expect_true(all(per_factor(part$freed$int) == 1L))
  expect_match(part$summary, "not supported")
})

# Item Order: three synthetic structures whose truth is known, and awkward input.
sim <- function(n = 1500, seed = 1) { set.seed(seed); n }
flip <- function(x, p = 0.01) as.integer(xor(x, runif(length(x)) < p))

test_that("a Guttman scale comes out as one line", {
  n <- sim()
  th <- rnorm(n)
  d <- as.data.frame(sapply(c(-1.2, -0.5, 0.2, 0.9), function(b) flip(th > b)))
  names(d) <- paste0("g", 1:4)
  x <- orderWiz(d, names(d))
  expect_equal(x$items$item, paste0("g", 1:4))
  expect_equal(x$items$level, 1:4)
  expect_match(x$reading[1], "one line")
  expect_equal(x$implications$from, c("g2", "g3", "g4"))
  expect_equal(x$implications$to, c("g1", "g2", "g3"))
  expect_equal(nrow(x$scales), 1)
  expect_equal(x$info$unscalable, character())
  expect_gt(x$fit$fitting[1], 0.9)
  expect_s3_class(plot(x), "ggplot")
})

test_that("a branch, a shared rung and an unrelated item are told apart", {
  n <- sim(seed = 2)
  s <- rnorm(n); e <- rnorm(n)
  base <- s > -0.3
  d <- data.frame(a1 = flip(base & runif(n) < .98), a2 = flip(base & runif(n) < .98),   # same rung
                  b = flip(base & s > 0.5),                                             # builds on a
                  c = flip(base & e > 0.6),                                             # builds on a, not on b
                  z = as.integer(runif(n) < 0.3))                                       # unrelated
  x <- orderWiz(d, names(d))
  it <- x$items
  expect_equal(it$level[match(c("a1", "a2", "b", "c"), it$item)], c(1, 1, 2, 2))
  expect_true(is.na(it$level[it$item == "z"]))
  expect_match(it$role[it$item == "z"], "apart")
  expect_equal(it$same.as[it$item == "a1"], "a2")
  expect_setequal(x$implications$from, c("b", "c"))
  expect_true(all(grepl("a1 = a2|a2 = a1", x$implications$to)))
  expect_match(x$reading[1], "branches")
  expect_true("z" %in% x$info$unscalable)
  # the order fits many more respondents than a line through all five items
  expect_gt(x$fit$fitting[2], x$fit$fitting[1] + 0.2)
  expect_output(print(x), "Apart from the order: z")
})

test_that("unrelated items give no order, and awkward input is handled", {
  n <- sim(seed = 3)
  d <- as.data.frame(matrix(as.integer(runif(n * 4) < 0.5), n))
  x <- orderWiz(d, names(d))
  expect_null(x$implications)
  expect_null(x$scales)
  expect_match(x$reading[1], "do not form an order")
  expect_s3_class(plot(x), "ggplot")

  d$const <- 1L
  d$likert <- sample(1:5, n, replace = TRUE)
  y <- orderWiz(d, names(d))
  expect_equal(y$info$k, 5)
  expect_true(any(grepl("same answer: const", y$notes)))
  expect_true(any(grepl("split at the midpoint", y$notes)))
  expect_error(orderWiz(d, c("V1", "V2")), "at least three")
  expect_error(orderWiz(d, c("V1", "const", "V2")), "Fewer than three items vary")
})

test_that("the example debt items give the structure they were built with, also in jamovi", {
  data(wizlab, package = "psyWiz")
  debt <- grep("^debt_", names(wizlab), value = TRUE)
  x <- orderWiz(wizlab, debt)
  lev <- stats::setNames(x$items$level, x$items$item)
  expect_equal(unname(lev[c("debt_family", "debt_friends", "debt_bank", "debt_lender", "debt_card")]), c(1, 1, 2, 2, 3))
  expect_true(is.na(lev[["debt_mortgage"]]))
  # (the jamovi functions read a bare name as a column name, so the names are spelled out)
  r <- itemOrder(wizlab, vars = paste0("debt_", c("family", "friends", "bank", "card", "lender", "mortgage")), showPatterns = TRUE, patternsN = 5)
  expect_equal(nrow(r$itemTable$asDF), 6)
  expect_equal(nrow(r$implications$asDF), 3)
  expect_equal(nrow(r$patterns$asDF), 5)
  expect_equal(r$hasse$state$kind, "hasse")
  expect_match(r$reading$content, "branches")
  expect_match(itemOrder(wizlab, vars = c("debt_family", "debt_friends"))$todo$content, "at least three")
  k <- orderWiz(wizlab, paste0("know", 1:6))
  expect_equal(k$items$level, 1:6)
})

context("Testing IE functionality")

test_that("Test IE function", {
  y <- c(rep(1, 10), rep(0, 10))
  phat <- seq(1, 0, length = length(y))
  Metric <- IE(y, phat)
  expect_equal(Metric, 2)
})
context("Testing RKL functionality")

test_that("Test RKL function", {
  y <- c(rep(1, 10), rep(0, 10))
  phat <- seq(1, 0, length = length(y))
  Metric <- RKL(y, phat)
  expect_equal(Metric, 10)
})
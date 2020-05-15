context("Testing AHR functionality")

test_that("Test AHR function", {
  y <- c(rep(1, 10), rep(0, 10))
  phat <- seq(1, 0, length = length(y))
  Metric <- AHR(y, phat)
  expect_equal(Metric, 1)
})
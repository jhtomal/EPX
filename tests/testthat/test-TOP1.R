context("Testing TOP1 functionality")

test_that("Test TOP1 function", {
  y <- c(rep(1, 10), rep(0, 10))
  phat <- seq(1, 0, length = length(y))
  Metric <- TOP1(y, phat)
  expect_equal(Metric, 1)
})
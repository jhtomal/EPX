context("Testing BNsample data")

test_that("Test BNsample data", {
  expect_equal(dim(BNsample), c(1000, 25))
})

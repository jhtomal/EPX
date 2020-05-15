context("Testing BNhold data")

test_that("Test BNhold data", {
  expect_equal(dim(BNhold), c(3946, 25))
})

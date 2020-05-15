context("Testing Harvest data")

test_that("Test Harvest data", {
  expect_equal(dim(harvest), c(190, 4))
})

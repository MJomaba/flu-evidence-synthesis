context("Cholesky")

test_that("Cholesky returns the lower decomposition matrix", {
  A <- matrix(c(4,12,-16,12,37,-43,-16,-43,98), ncol = 3)
  B <- cholesky(A)
  expect_identical( B, matrix(c(2,6,-8,0,1,5,0,0,3),ncol = 3) )
})

test_that("bounds_ate argument checks work", {
  x <- data.frame(X = 1, A = 1, Y = c(1, 0))
  x_nonbinary <- data.frame(X = 1, A = 1, Y = c(0.5, 0))

  expect_error(ate_bounds())
  expect_error(ate_bounds(data = 1))
  expect_error(ate_bounds(data = x))
  expect_error(ate_bounds(data = x, X = 1, A = 2, Y = 3))
  expect_error(ate_bounds(data = x, c("X1", "X2"), "A", "Y"))
  expect_error(ate_bounds(data = x, "X", "A1", "Y"))
  expect_error(ate_bounds(data = x, "X", "A", "Y1"))
  expect_error(ate_bounds(data = x_nonbinary, "X", "A", "Y"))
  expect_error(ate_bounds(data = x, "X", "A", "Y", inner_folds = 0))
})

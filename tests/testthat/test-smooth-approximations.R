test_that("smooth approximation of I[x > t] passes spot checks", {
  expect_equal(s_gt(1, 0, 0), 1)  # s_g(1, 0, 0) = I[1 > 0]
  expect_equal(s_gt(-1, 0, 0), 0)  # s_g(-1, 0, 0) = I[-1 > 0]

  expect_equal(s_gt(0, 0, 1), 0)
  expect_equal(s_gt(1, 0, 1), 1)
  expect_lt(s_gt(0.5, 0, 1), 1)
  expect_gt(s_gt(0.5, 0, 1), 0)
})

test_that("smooth approximation of I[x < t] passes spot checks", {
  expect_equal(s_lt(1, 0, 0), 0)  # s_l(1, 0, 0) = I[1 < 0]
  expect_equal(s_lt(-1, 0, 0), 1)  # s_l(-1, 0, 0) = I[-1 < 0]

  expect_equal(s_lt(0, 0, 1), 0)
  expect_equal(s_lt(-1, 0, 1), 1)
  expect_lt(s_lt(-0.5, 0, 1), 1)
  expect_gt(s_lt(-0.5, 0, 1), 0)
})

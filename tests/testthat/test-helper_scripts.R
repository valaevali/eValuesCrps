library(testthat)

test_that("test forecast_input", {
  expect_type(
    forecast_input("perfect", rnorm(100), 1),
    "list"
  )
})

test_that("test forecast_input mixnorm", {
  expect_type(
    forecast_input("mixnorm", "mu" = cbind(rnorm(100), rnorm(100) + 1), "sd" = matrix(nrow = 100, ncol = 2, 1), "w" = matrix(nrow = 100, ncol = 2, 1 / 2)),
    "list"
  )
})

test_that("test optim_inf_values returns inf_value", {
  expect_equal(
    optim_inf_value(\(x){x^2}, min.value=-2, max.value=2),
    0
  )
})

library(testthat)

test_that("test forecast_input", {
  expect_type(
    forecast_input(rnorm(100), 1),
    "list"
  )
})

test_that("test forecast_input mixnorm", {
  expect_type(
    forecast_input("mu" = cbind(rnorm(100), stats::rnorm(100) + 1), "sd" = matrix(nrow = 100, ncol = 2, 1), "w" = matrix(nrow = 100, ncol = 2, 1 / 2)),
    "list"
  )
})

test_that("test optim_inf_values returns inf_value", {
  expect_equal(
    optim_inf_value(\(x){x^2}, min.value=-2, max.value=2),
    0
  )
})

mu <- stats::rnorm(10)
inf.fun <- \(x, j) {scoringRules::crps_norm(y = x, mean = mu[j], sd = 1)}
test_that("input for perfect returns correct inf fun", {
  expect_equal(
    create_crps_fun(mu = mu, sd = 1)$inf.fun(0.5, 5),
    inf.fun(0.5, 5)
  )
})

inf.fun <- \(x, j) {scoringRules::crps_norm(y = x, mean = 0, sd = 2)}
test_that("input for climatological returns correct inf fun", {
  expect_equal(
    create_crps_fun(mu = 0, sd = 2)$inf.fun(0.5, 5),
    inf.fun(0.5, 5)
  )
})

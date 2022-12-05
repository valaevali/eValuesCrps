library(testthat)

test_that("test e_value does return evalue", {
  expect_type(
    e_value(rnorm(10), forecast_input(mu = rnorm(10), sd = 1), forecast_input(mu = -rnorm(10), sd = 1)),
    "list"
  )
})

mu <- rnorm(10)
tau <- sample(c(-1, 1), 10, replace = TRUE)
test_that("test e_value does return evalue for mixnorm", {
  expect_type(
    e_value(rnorm(10), forecast_input(mu = mu, sd = 1), forecast_input(mu = cbind(mu, mu + tau), sd = matrix(nrow = 10, ncol = 2, 1), w = matrix(nrow = 10, ncol = 2, 1 / 2))),
    "list"
  )
})

test_that("test forecast_input", {
  expect_type(
    forecast_input(rnorm(100), 1),
    "list"
  )
})

test_that("test forecast_input mixnorm", {
  expect_type(
    forecast_input("mu" = cbind(rnorm(100), rnorm(100) + 1), "sd" = matrix(nrow = 100, ncol = 2, 1), "w" = matrix(nrow = 100, ncol = 2, 1 / 2)),
    "list"
  )
})

test_that("test optim_inf_values returns inf_value", {
  expect_equal(
    optim_inf_value(\(x){x^2}, min.value=-2, max.value=2),
    0
  )
})

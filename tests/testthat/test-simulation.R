library(testthat)

test_that("test crps_sim for e_value", {
  expect_s3_class(
    sim_e_values(n.it = 2, n.obs = 2),
    "list"
  )
})
#' This method provides predefined forecasts.
#' @export
special_forecasts <- function(mu, tau, n.obs, loosing.power.forecasts = FALSE, usual.forecasts = TRUE) {
  forecasts <- list("perfect" = list("mu" = mu, "sd" = 1, "main" = TRUE))
  if (loosing.power.forecasts == TRUE) {
    forecasts <- base::append(forecasts, list(
      "perfect-m-0.01" = list("mu" = mu + 0.01, "sd" = 1),
      "perfect-m-0.025" = list("mu" = mu + 0.025, "sd" = 1),
      "perfect-m-0.05" = list("mu" = mu + 0.05, "sd" = 1),
      "perfect-m-0.1" = list("mu" = mu + 0.1, "sd" = 1),
      "perfect-m-0.2" = list("mu" = mu + 0.2, "sd" = 1),
      "perfect-m-0.3" = list("mu" = mu + 0.3, "sd" = 1),
      "perfect-m-0.4" = list("mu" = mu + 0.4, "sd" = 1),
      "perfect-m-0.5" = list("mu" = mu + 0.5, "sd" = 1),
      "perfect-m-0.6" = list("mu" = mu + 0.6, "sd" = 1),
      "perfect-m-0.7" = list("mu" = mu + 0.7, "sd" = 1),
      "perfect-m-0.8" = list("mu" = mu + 0.8, "sd" = 1),
      "perfect-m-0.9" = list("mu" = mu + 0.9, "sd" = 1),
      "perfect-m-1.0" = list("mu" = mu + 1.0, "sd" = 1),
      "perfect-s-0.01" = list("mu" = mu, "sd" = 1 + 0.01),
      "perfect-s-0.025" = list("mu" = mu, "sd" = 1 + 0.025),
      "perfect-s-0.05" = list("mu" = mu, "sd" = 1 + 0.05),
      "perfect-s-0.1" = list("mu" = mu, "sd" = 1 + 0.1),
      "perfect-s-0.2" = list("mu" = mu, "sd" = 1 + 0.2),
      "perfect-s-0.3" = list("mu" = mu, "sd" = 1 + 0.3),
      "perfect-s-0.4" = list("mu" = mu, "sd" = 1 + 0.4),
      "perfect-s-0.5" = list("mu" = mu, "sd" = 1 + 0.5),
      "perfect-s-0.6" = list("mu" = mu, "sd" = 1 + 0.6),
      "perfect-s-0.7" = list("mu" = mu, "sd" = 1 + 0.7),
      "perfect-s-0.8" = list("mu" = mu, "sd" = 1 + 0.8),
      "perfect-s-0.9" = list("mu" = mu, "sd" = 1 + 0.9),
      "perfect-s-1.0" = list("mu" = mu, "sd" = 1 + 1.0)
    ))
  }

  if (usual.forecasts == TRUE) {
    forecasts <- base::append(forecasts, list(
      "climatological" = list("mu" = 0, "sd" = 2),
      "sign-reversed" = list("mu" = -mu, "sd" = 1),
      "unfocused" = list("mu" = cbind(mu, mu + tau), "sd" = matrix(nrow = n.obs, ncol = 2, 1), "w" = matrix(nrow = n.obs, ncol = 2, 1 / 2))
    ))
  }
  forecasts
}


#' This method provides predefined forecasts for the log score.
#' @export
special_forecasts_logs <- function(mu, tau, n.obs, loosing.power.forecasts = FALSE, usual.forecasts = TRUE) {
  forecasts <- list("perfect" = list("mu" = mu, "sd" = 1, "main" = TRUE, "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1), nrow = n.obs)}, "inf.fun" = \(x,j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1)}))
  if (loosing.power.forecasts == TRUE) {
    forecasts <- base::append(forecasts, list(
      "perfect-m-0.01" =  list("mu" = mu + 0.01,  "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.01,  sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.01,  sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.01, sd = 1)} ),
      "perfect-m-0.025" = list("mu" = mu + 0.025, "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.025, sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.025, sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.025,sd = 1)} ),
      "perfect-m-0.05" =  list("mu" = mu + 0.05,  "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.05,  sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.05,  sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.05, sd = 1)} ),
      "perfect-m-0.1" =   list("mu" = mu + 0.1,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.1,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.1,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.1,  sd = 1)} ),
      "perfect-m-0.2" =   list("mu" = mu + 0.2,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.2,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.2,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.2,  sd = 1)} ),
      "perfect-m-0.3" =   list("mu" = mu + 0.3,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.3,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.3,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.3,  sd = 1)} ),
      "perfect-m-0.4" =   list("mu" = mu + 0.4,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.4,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.4,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.4,  sd = 1)} ),
      "perfect-m-0.5" =   list("mu" = mu + 0.5,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.5,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.5,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.5,  sd = 1)} ),
      "perfect-m-0.6" =   list("mu" = mu + 0.6,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.6,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.6,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.6,  sd = 1)} ),
      "perfect-m-0.7" =   list("mu" = mu + 0.7,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.7,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.7,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.7,  sd = 1)} ),
      "perfect-m-0.8" =   list("mu" = mu + 0.8,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.8,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.8,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.8,  sd = 1)} ),
      "perfect-m-0.9" =   list("mu" = mu + 0.9,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 0.9,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 0.9,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 0.9,  sd = 1)} ),
      "perfect-m-1.0" =   list("mu" = mu + 1.0,   "sd" = 1,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean =  mu + 1.0,   sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu + 1.0,   sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean =  mu[j] + 1.0,  sd = 1)} ),
      "perfect-s-0.01" =  list("mu" = mu, "sd" = 1 + 0.01, "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.01  )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.01 ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.01 )} ),
      "perfect-s-0.025" = list("mu" = mu, "sd" = 1 + 0.025,"method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.025 )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.025), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.025)} ),
      "perfect-s-0.05" =  list("mu" = mu, "sd" = 1 + 0.05, "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.05  )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.05 ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.05 )} ),
      "perfect-s-0.1" =   list("mu" = mu, "sd" = 1 + 0.1,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.1   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.1  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.1  )} ),
      "perfect-s-0.2" =   list("mu" = mu, "sd" = 1 + 0.2,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.2   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.2  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.2  )} ),
      "perfect-s-0.3" =   list("mu" = mu, "sd" = 1 + 0.3,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.3   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.3  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.3  )} ),
      "perfect-s-0.4" =   list("mu" = mu, "sd" = 1 + 0.4,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.4   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.4  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.4  )} ),
      "perfect-s-0.5" =   list("mu" = mu, "sd" = 1 + 0.5,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.5   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.5  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.5  )} ),
      "perfect-s-0.6" =   list("mu" = mu, "sd" = 1 + 0.6,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.6   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.6  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.6  )} ),
      "perfect-s-0.7" =   list("mu" = mu, "sd" = 1 + 0.7,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.7   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.7  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.7  )} ),
      "perfect-s-0.8" =   list("mu" = mu, "sd" = 1 + 0.8,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.8   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.8  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.8  )} ),
      "perfect-s-0.9" =   list("mu" = mu, "sd" = 1 + 0.9,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 0.9   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 0.9  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 0.9  )} ),
      "perfect-s-1.0" =   list("mu" = mu, "sd" = 1 + 1.0,  "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = mu, sd = 1 + 1.0   )}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = mu, sd = 1 + 1.0  ), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = mu[j], sd = 1 + 1.0  )} )
    ))
  }

  if (usual.forecasts == TRUE) {
    forecasts <- base::append(forecasts, list(
      "climatological" = list("mu" = 0, "sd" = 2, "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = 0, sd = 2)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = 0, sd = 2), nrow = n.obs)}),
      "sign-reversed" = list("mu" = -mu, "sd" = 1, "method" = 'norm-logs', "fun" = \(y) {scoringRules::logs_norm(y = y, mean = -mu, sd = 1)}, "sample.fun" = \(n, n.obs) {matrix(rnorm(n * n.obs, mean = -mu, sd = 1), nrow = n.obs)}, "inf.fun" = \(x, j) {scoringRules::logs_norm(y = x, mean = -mu[j], sd = 1)} ),
      "unfocused" = list("mu" = cbind(mu, mu + tau), "sd" = matrix(nrow = n.obs, ncol = 2, 1), "w" = matrix(nrow = n.obs, ncol = 2, 1 / 2), "method" = 'mixnorm-logs',
                         "fun" = \(y) { scoringRules::logs_mixnorm(y = y, m = cbind(mu, mu + tau), s = matrix(nrow = n.obs, ncol = 2, 1), w = matrix(nrow = n.obs, ncol = 2, 1 / 2)) },
                         "crps.fun.y.matrix" = \(y) { sapply(1:dim(y)[2], \(i) {scoringRules::logs_mixnorm(y = y[, i], m = cbind(mu, mu + tau), s = matrix(nrow = n.obs, ncol = 2, 1), w = matrix(nrow = n.obs, ncol = 2, 1 / 2))}) },
                         "sample.fun" = \(n, n.obs) { matrix(stats::rnorm(n * n.obs, mean = cbind(mu, mu + tau) * matrix(nrow = n.obs, ncol = 2, 1 / 2), sd = matrix(nrow = n.obs, ncol = 2, 1) * matrix(nrow = n.obs, ncol = 2, 1 / 2)^2), nrow = n.obs) },
                         "inf.crps.fun" = \(x, j) { scoringRules::logs_mixnorm(y = x, m = as.matrix(t(cbind(mu, mu + tau)[j,]), nrow = 1), s = as.matrix(t(matrix(nrow = n.obs, ncol = 2, 1)[1,])), w = as.matrix(t(matrix(nrow = n.obs, ncol = 2, 1 / 2)), nrow = n.obs)[1,]) } )
    ))
  }
  forecasts
}
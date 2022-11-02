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


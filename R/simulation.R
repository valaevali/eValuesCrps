#' This method performs an simulation for e-values and their power.
#'
#' @param n.it = 200, number of iterations the simulation should make.
#' @param n.obs = 100, number of simulated observation in each iteration step.
#' @param loosing.power.only = FALSE, if "TRUE" only the e-values for comparing forecasts with a bias in the mean or variance are calculated.
#' @param method = list("GRAPA", "lambda", "alternative", "alternative-mean"), is a list containing all the method names for calculating the different lambdas.
#' @param lambda = 0.5, lambda entry for a fixed value.
#' @param p.value.method = "t", can be "t" for t-test and "dm" for dm.test of the package forecast.
#' @param forecasts.input = NA, one can pass a list for different forecasts to be taken into account for calculating the e-values, they have to be of the form list("additional-forecast" = list("mu" = mu, "sd" = 1)) or for mixed forecasts use "unfocused" = list("mu" = cbind(mu, mu + tau), "sd" = matrix(nrow = n.obs, ncol = 2, 1), "w" = matrix(nrow = n.obs, ncol = 2, 1 / 2)). Use [forecast_input()] to produce correct form.
#' @param loosing.power.forecasts = loosing.power.only, set to TRUE if you want to include the predefined forecasts with a bias in the mean or variance.
#' @param usual.forecasts = TRUE, set to FALSE if you want to exclude the predefined usual forecasts.
#' @param file.folder = getwd(), where to save the output to.

#' @export
sim_e_values <- function(n.it = 200, n.obs = 100, loosing.power.only = FALSE,
                         method = list("GRAPA", "lambda", "alternative"), lambda = 0.5, p.value.method = "t",
                         forecasts.input = NA, loosing.power.forecasts = loosing.power.only, usual.forecasts = TRUE, file.folder = getwd()) {
  checkedInput <- check_input_simulation(n.it, n.obs, loosing.power.only, method, lambda, p.value.method, forecasts.input, loosing.power.forecasts, usual.forecasts, file.folder)
  n.it <- checkedInput$n.it
  n.obs <- checkedInput$n.obs
  loosing.power.only <- checkedInput$loosing.power.only
  method <- checkedInput$method
  lambda <- checkedInput$lambda
  p.value.method <- checkedInput$p.value.method
  forecasts.input <- checkedInput$forecasts.input
  loosing.power.forecasts <- checkedInput$loosing.power.forecasts
  usual.forecasts <- checkedInput$usual.forecasts
  file.folder <- checkedInput$file.folder

  start.time <- Sys.time()
  cl <- parallel::makeCluster(parallel::detectCores() - 2, type = "PSOCK")
  doSNOW::registerDoSNOW(cl = cl)
  on.exit(parallel::stopCluster(cl = cl))
  pb <- progress::progress_bar$new(total = n.it + 1)

  # allowing progress bar to be used in dopar
  opts <- list(progress = \() { pb$tick(1) })

  pb$tick(1)
  result <- foreach::foreach(i = seq(1, n.it), .combine = \(x, y) { data.table::rbindlist(list(x, y)) }, .options.snow = opts, .packages = "dplyr") %dopar% {
    sapply(list.files(paste0(getwd(), "/R"), pattern = "*.R", full.names = TRUE, recursive = FALSE), \(f) { source(f) })

    mu <- rnorm(n.obs)
    y <- rnorm(n.obs, mu)
    tau <- sample(c(-1, 1), n.obs, replace = TRUE)

    forecasts <- base::append(forecasts.input, special_forecasts(mu, tau, n.obs, loosing.power.forecasts, usual.forecasts))
    forecasts <- forecasts[!is.na(forecasts)]

    temp <- outer(forecasts, forecasts, Vectorize(\(f, g) {
      if (!identical(f, g) && ifelse(loosing.power.only, ifelse(!is.null(g$main), TRUE, FALSE), TRUE)) {
        return(e_value(y = y, crps.F.para = f, crps.G.para = g, it = i, method = method, lambda = lambda, p.value.method = p.value.method))
      }
    }))
    temp.tibble <- tibble::as_tibble(temp)
    temp.tibble <- temp.tibble %>%
      mutate(names.F = rownames(temp)) %>%
      tidyr::pivot_longer(cols = colnames(temp.tibble), names_to = "names.G") %>%
      tidyr::unnest_wider(value) %>%
      na.omit(p.value)
    temp.tibble
  }

  timestamp <- format(Sys.time(), format = "%Y-%m-%dT%H-%M-%S")
  result.eval <- result %>%
    group_by(names.F, names.G) %>%
    mutate(across(contains(c("p.value", "prod")), ~100 / n.it * sum(.x <= 0.05), .names = "{.col}.H0.rej")) %>%
    select(contains(c("names", "H0"))) %>%
    distinct() %>%
    arrange(names.F, names.G)

  print(Sys.time() - start.time)
  result.fin <- list("uncompacted" = result, "evaluated" = result.eval)

  saveRDS(result.fin, paste0(file.folder, "/target/run-", n.obs, "-", n.it, "-", timestamp, ".rds"))
  if (loosing.power.only) {
    pdf(paste0(getwd(), "/target/plot-", n.obs, "-", n.it, "-", timestamp, ".pdf"))
    print_rej_rate_perfect_loosing_power(result.fin)
    dev.off()
  }

  return(result.fin)
}
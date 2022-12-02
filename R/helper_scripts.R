# Title     : Helper Scripts
# Created by: Velerie Haftka
# Created on: 26.01.2022

check_input_simulation <- function(n.it, n.obs, loosing.power.only, method, lambda, p.value.method, forecasts.input, loosing.power.forecasts, usual.forecasts, file.folder) {
  if (is.na(n.it) || !is.numeric(n.it)) {
    warning("Number of iterations should be a not-null numeric. Setting to default 200.")
    n.it <- 200
  }
  if (is.na(n.obs) || !is.numeric(n.obs)) {
    warning("Number of observations should be a not-null numeric. Setting to default 200.")
    n.obs <- 200
  }
  if (is.na(loosing.power.only) || !is.logical(loosing.power.only)) {
    warning("Parameter 'loosing.power.only' should be a not-null logical. Setting to default FALSE.")
    loosing.power.only <- FALSE
  }
  if (!("GRAPA" %in% method) && !("lambda" %in% method) && !("alternative" %in% method) && !("alternative-mean" %in% method)) {
    warning("Parameter 'method' did not match any defined methods, setting to default 'lambda' for a fixed value")
    method <- list("lambda")
  }
  if (is.na(lambda) && ("lambda" %in% method)) {
    warning("Method with fixed lambda is choosen, but no lambda is provided. Setting to default lambda = 0.5")
    lambda <- 0.5
  }
  if (is.na(p.value.method)) {
    warning("Parameter 'p.value.method' is na, hence no p.value will be calculated")
  }
  if (!is.na(p.value.method) &&
    !("t" %in% p.value.method) &&
    !("dm" %in% p.value.method)) {
    warning("Parameter 'p.value.method' has to be one of ('t','dm'). Setting to default 't'")
    p.value.method <- "t"
  }
  if (!is.na(forecasts.input) && !is.list(forecasts.input)) {
    warning("Parameter 'forecasts.input' should be a list of lists. Ignoring input.")
    forecasts.input <- NA
  }
  if (is.na(loosing.power.forecasts) || !is.logical(loosing.power.forecasts)) {
    warning("Parameter 'loosing.power.forecasts' should be a not-null logical. Setting to default FALSE.")
    loosing.power.forecasts <- FALSE
  }
  if (is.na(usual.forecasts) || !is.logical(usual.forecasts)) {
    warning("Parameter 'usual.forecasts' should be a not-null logical. Setting to default TRUE.")
    usual.forecasts <- TRUE
  }
  if (is.na(file.folder)) {
    warning("Parameter 'file.folde' should not be null. Setting to default: getwd().")
    file.folder <- getwd()
  }

  return(list("n.it" = n.it, "n.obs" = n.obs, "loosing.power.only" = loosing.power.only, "method" = method, "lambda" = lambda, "p.value.method" = p.value.method, "forecasts.input" = forecasts.input, "loosing.power.forecasts" = loosing.power.forecasts, "usual.forecasts" = usual.forecasts, "file.folder" = file.folder))
}

check_input <- function(y, crps.F.para, crps.G.para, it, method, lambda, p.value.method) {
  if (is.na(y) || !is.vector(y)) {
    stop("The parameter y had to be a vector.")
  }

  check_input_crps_para(crps.F.para, "F")
  check_input_crps_para(crps.G.para, "G")

  if (is.na(it)) {
    warning("Setting iteration to '1'")
    it <- 1
  }
  if (!("GRAPA" %in% method) && !("lambda" %in% method) && !("alternative" %in% method) && !("alternative-mean" %in% method)) {
    warning("Parameter 'method' did not match any defined methods, setting to default 'lambda' for a fixed value")
    method <- list("lambda")
  }

  if (is.na(lambda) && ("lambda" %in% method)) {
    warning("Method with fixed lambda is choosen, but no lambda is provided. Setting to default lambda = 0.5")
    lambda <- 0.5
  }

  if (is.na(p.value.method)) {
    warning("Parameter 'p.value.method' is na, hence no p.value will be calculated")
  }

  if (!is.na(p.value.method) &&
    !("t" %in% p.value.method) &&
    !("dm" %in% p.value.method)) {
    warning("Parameter 'p.value.method' has to be one of ('t','dm'). Setting to default 't'")
    p.value.method <- "t"
  }

  return(list("it" = it, "method" = method, "lambda" = lambda, "p.value.method" = p.value.method))
}

check_input_crps_para <- function(crps.para, f.g) {
  if (all(is.na(crps.para))) {
    stop(sprintf("The CRPS parameter should not be empty for %s.", f.g))
  }
  if (!("mu" %in% names(crps.para))) {
    stop(sprintf("The CRPS parameter 'mu' should not be empty for %s.", f.g))
  }
  if (!("sd" %in% names(crps.para))) {
    stop(sprintf("The CRPS parameter 'sd' should not be empty for %s.", f.g))
  }
  if (!is.numeric(crps.para$mu) && !is.vector(crps.para$mu)) {
    stop(sprintf("The CRPS parameter 'mu' should be a constant or a vector for %s.", f.g))
  }
  if (!is.numeric(crps.para$sd) && !is.vector(crps.para$sd)) {
    stop(sprintf("The CRPS parameter 'sd' should be a constant or a vector for %s.", f.g))
  }
  if ("w" %in% names(crps.para) &&
    !is.na(crps.para$w) &&
    !is.numeric(crps.para$w) &&
    !(is.vector(crps.para$w) && !is.matrix(crps.para$w))) {
    stop(sprintf("The CRPS parameter 'w' should be either NA or a constant or a vector or a matrix for %s.", f.g))
  }

  if ((is.matrix(crps.para$mu) && (!is.matrix(crps.para$sd) || ("w" %in% names(crps.para) && !is.matrix(crps.para$w)))) ||
    (is.matrix(crps.para$sd) && (!is.matrix(crps.para$mu) || ("w" %in% names(crps.para) && !is.matrix(crps.para$w)))) ||
    (("w" %in% names(crps.para) && !is.matrix(crps.para$w)) && (!is.matrix(crps.para$sd) || !is.matrix(crps.para$mu)))) {
    stop(sprintf("One of the parameter for '%s' is a matrix, but not the others. Make sure all the parameters for one forecast have the same form", f.g))
  }
}

#' This is a helper function which returns the correct form for the forecast input for [sim_e_values()].
#'
#' @param forecast.name is the name of the forecast.
#' @param mu is the mean of the forecast, can be constant, vector or matrix.
#' @param sd is the variance of the forecast, can be a constant, vector or matrix.
#' @param w = NA, optional, is the weight of a mixed forecast, can be a constant, vector or matrix.
#' @export
forecast_input <- function(forecast.name, mu, sd, w = NA) {
  if (is.na(w)) {
    return(list(forecast.name = list("mu" = mu, "sd" = sd)))
  }
  return(list(forecast.name = list("mu" = mu, "sd" = sd, "w" = w)))
}

# This method makes a gitter search for global minima
optim_inf_value <- function(f, start.points = 2, min.value = 0.0001, max.value = 1) {
  min(sapply(seq(min.value, max.value, length.out = start.points), \(y) { optim(y, f, lower = min.value, upper = max.value, method = "L-BFGS-B")$value }))
}

#' This method returns the R object save it the file with getwd() + filename.
#' @param filename is the name of the file.
#' @returns the R object in the file.
#' @examples \dontrun{getFile("/test.rds")}
#' @export
getFile <- function(filename) {
  readRDS(paste0(getwd(), filename))
}

#' This methods creates the additional functions from the input parameters for the calculations of lambda.
#' @param n.obs = 200, number of observations.
#' @param mu = 0, the mean of the forecast.
#' @param sd = 1, the variance of the forecast.
#' @param w = 1, the weight of a mixed forecast.
#' @param ... Additional parameters which are ignored.
#' @returns A list of method = ("norm"|"mixnorm"), fun = \(y) {scoringRules::crps_norm(y=y, mean=mu, sd=sd} or \(y) { scoringRules::crps_mixnorm(y = y, m = mu, s = sd, w = w) } for mixed norm,
#'          crps.fun.y.matrix is a function only special if method = "mixnorm" to calculate the alternative betting adaptive to forecasts otherwise it is equal to fun,
#'          rnorm is the function to calculate randomly normally distributed variables, inf.fun is the function to calculate the infimum for this forecast.
#' @examples \dontrun{create_crps_fun(n.obs = 200, mu = rnorm(200), sd = 1)}
#' @export
create_crps_fun <- function(n.obs = 200, mu = 0, sd = 1, w = 1, ...) {
  if (is.matrix(mu) || is.matrix(sd) || is.matrix(w)) {
    method <- 'mixnorm'
    crps.fun <- \(y) { scoringRules::crps_mixnorm(y = y, m = mu, s = sd, w = w) }
    crps.fun.y.matrix <- \(y) { sapply(1:dim(y)[2], \(i) { scoringRules::crps_mixnorm(y = y[, i], m = mu, s = sd, w = w) }) }
    rnorm.fun <- \(n) { matrix(rnorm(n * n.obs, mean = mu, sd = sd), nrow = n.obs) }
    inf.crps.fun <- \(x, j) { scoringRules::crps_mixnorm(y = x, m = as.matrix(t(mu[j,]), nrow = 1), s = as.matrix(t(sd[1, ])), w = as.matrix(t(w[1, ]))) }
  } else {
    method <- 'norm'
    crps.fun <- \(y) { scoringRules::crps_norm(y = y, mean = mu, sd = sd) }
    crps.fun.y.matrix <- crps.fun
    rnorm.fun <- \(n) { matrix(rnorm(n * n.obs, mean = mu, sd = sd), nrow = n.obs) }
    inf.crps.fun <- \(x, j) { scoringRules::crps_norm(y = x, mean = mu[j], sd = sd) }
  }
  return(list("method" = method, "fun" = crps.fun, "crps.fun.y.matrix" = crps.fun.y.matrix,
              "rnorm" = rnorm.fun, "inf.fun" = inf.crps.fun))
}

#' This method creates a pdf with the plots of the rejection rate of the forecasts with bias in the mean or variance.
#' @param dt is the result of the [sim_e_values(loosing.power.forecasts = TRUE)]. But it has to have the loosing.power.forecasts = TRUE.
#' @export
print_rej_rate_perfect_loosing_power <- function(dt) {
  to.print <- dt$evaluated %>%
    filter(grepl("perfect", names.F) & names.G == 'perfect') %>%
    mutate(e = stringr::str_extract(names.F, "[.0-9]+"),
           mean.sd = stringr::str_extract(names.F, "perfect-[a-z]")
    ) %>%
    ungroup() %>%
    select(-c(names.F, names.G)) %>%
    tidyr::pivot_longer(!c(e, mean.sd), names_to = "key", values_to = "rej_rate") %>%
    mutate(key = stringr::str_replace_all(key, ".prod.H0.rej", "")) %>%
    mutate(key = stringr::str_replace_all(key, ".H0.rej", "")) %>%
    arrange(key)

  g.m <- ggplot2::ggplot(to.print %>% filter(grepl("-m", mean.sd)) %>% select(-mean.sd), ggplot2::aes(x = as.numeric(e), y = rej_rate, color = key)) +
    ggplot2::geom_line(ggplot2::aes(group = key)) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dotted") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    ggplot2::labs(x = "epsilon") +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.s <- ggplot2::ggplot(to.print %>% filter(grepl("-s", mean.sd)) %>% select(-mean.sd), ggplot2::aes(x = as.numeric(e), y = rej_rate, color = key)) +
    ggplot2::geom_line(ggplot2::aes(group = key)) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dotted") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    ggplot2::labs(x = "epsilon") +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  Rmisc::multiplot(g.m, g.s, cols = 1)
}

#' This method creates a pdf with the plots of the histograms of the overall time difference of the e-values.
#' @param dt is the result of the [sim_e_values(loosing.power.forecasts = TRUE)]. But it has to have the loosing.power.forecasts = TRUE.
#' @export
print_e_values_histogram_loosing_power <- function(dt) {
  t.e.values <- dt$uncompacted %>%
    filter(grepl("perfect", names.F) & names.G == 'perfect') %>%
    select(contains(c("names", "e.value"))) %>%
    select(!contains(".prod")) %>%
    tidyr::unnest(contains("e.value")) %>%
    mutate(names = paste0(names.F, ".", names.G)) %>%
    select(-c(names.F, names.G)) %>%
    group_by(names) %>%
    dplyr::mutate(across(contains("e.value"), ~mean(.x), .names = "{.col}.mean"))

  g.grapa.m <- ggplot2::ggplot((t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(x = e.value.grapa, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(xintercept = e.value.grapa.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.grapa.s <- ggplot2::ggplot((t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(x = e.value.grapa, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(xintercept = e.value.grapa.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.05.m <- ggplot2::ggplot((t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(x = e.value.lambda, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(xintercept = e.value.lambda.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.05.s <- ggplot2::ggplot((t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(x = e.value.lambda, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(xintercept = e.value.lambda.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.cons.m <- ggplot2::ggplot((t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(x = e.value.alt.conf, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(xintercept = e.value.alt.conf.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.cons.s <- ggplot2::ggplot((t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(x = e.value.alt.conf, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(xintercept = e.value.alt.conf.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.conf.m <- ggplot2::ggplot((t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(x = e.value.alt.cons, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-m-", names))), ggplot2::aes(xintercept = e.value.alt.cons.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.conf.s <- ggplot2::ggplot((t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(x = e.value.alt.cons, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.e.values %>% filter(grepl("-s-", names))), ggplot2::aes(xintercept = e.value.alt.cons.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  Rmisc::multiplot(g.05.m, g.05.s, cols = 1)
  Rmisc::multiplot(g.grapa.m, g.grapa.s, cols = 1)
  Rmisc::multiplot(g.alt.conf.m, g.alt.conf.s, cols = 1)
  Rmisc::multiplot(g.alt.cons.m, g.alt.cons.s, cols = 1)
}

#' This method creates a pdf with the plots of the histograms of the overall time difference of the e-values.
#' @param dt is the result of the [sim_e_values(usual.forecasts = TRUE)]. But it has to have the usual.forecasts = TRUE.
#' @export
print_e_values_histogram_usual_forecasts <- function(dt) {
  t.e.values <- filter_tibble_for_usual_forecasts(dt$uncompacted) %>%
    select(contains(c("names", "e.value"))) %>%
    select(!contains(".prod")) %>%
    tidyr::unnest(contains("e.value")) %>%
    mutate(names = paste0(names.F, ".", names.G)) %>%
    select(-c(names.F, names.G)) %>%
    group_by(names) %>%
    dplyr::mutate(across(contains("e.value"), ~mean(.x), .names = "{.col}.mean"))

  g.grapa.m <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.grapa, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.grapa.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.grapa.s <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.grapa, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.grapa.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.05.m <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.lambda, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.lambda.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.05.s <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.lambda, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.lambda.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.cons.m <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.alt.conf, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.alt.conf.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.cons.s <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.alt.conf, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.alt.conf.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.conf.m <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.alt.cons, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.alt.cons.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.alt.conf.s <- ggplot2::ggplot(t.e.values, ggplot2::aes(x = e.value.alt.cons, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.e.values, ggplot2::aes(xintercept = e.value.alt.cons.mean, color = names)) +
    ggplot2::xlim(-1, 3) +
    ggplot2::ylim(0, 100000) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  Rmisc::multiplot(g.05.m, g.05.s, cols = 1)
  Rmisc::multiplot(g.grapa.m, g.grapa.s, cols = 1)
  Rmisc::multiplot(g.alt.conf.m, g.alt.conf.s, cols = 1)
  Rmisc::multiplot(g.alt.cons.m, g.alt.cons.s, cols = 1)
}

#' This method creates a pdf with the plots of the histograms of the difference of the crps values.
#' @param dt is the result of the [sim_e_values(loosing.power.forecasts = TRUE)]. But it has to have the loosing.power.forecasts = TRUE.
#' @export
print_crps_diff_histogram_loosing_power <- function(dt) {
  t.crps.dif <- dt$uncompacted %>%
    filter(grepl("perfect", names.F) & names.G == 'perfect') %>%
    select(names.F, crps.F, crps.G) %>%
    tidyr::unnest(c(crps.F, crps.G)) %>%
    mutate(dif.crps = crps.F - crps.G) %>%
    select(names.F, dif.crps) %>%
    group_by(names.F) %>%
    dplyr::mutate(mean = mean(dif.crps)) %>%
    ungroup()

  g.m <- ggplot2::ggplot(t.crps.dif %>% filter(grepl("-m-", names.F)), ggplot2::aes(x = dif.crps, fill = names.F, color = names.F)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.crps.dif %>% filter(grepl("-m-", names.F))), ggplot2::aes(xintercept = mean, color = names.F)) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  g.s <- ggplot2::ggplot(t.crps.dif %>% filter(grepl("-s-", names.F)), ggplot2::aes(x = dif.crps, fill = names.F, color = names.F)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = (t.crps.dif %>% filter(grepl("-s-", names.F))), ggplot2::aes(xintercept = mean, color = names.F)) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  Rmisc::multiplot(g.m, g.s, cols = 1)
}

#' This method creates a pdf with the plots of the histograms of the difference of the crps values.
#' @param dt is the result of the [sim_e_values(usual.forecasts = TRUE)]. But it has to have the usual.forecasts = TRUE.
#' @param filename to save the pdf to.
#' @examples \dontrun{print_crps_diff_histogram_usual_forecasts(dt, "/test")}
#' @export
print_crps_diff_histogram_usual_forecasts <- function(dt, filename) {
  t.crps.dif <- filter_tibble_for_usual_forecasts(dt$uncompacted) %>%
    select(names.F, names.G, it, crps.F, crps.G) %>%
    tidyr::unnest(c(crps.F, crps.G)) %>%
    mutate(dif.crps = crps.F - crps.G, names = paste0(names.F, "-", names.G)) %>%
    select(names, dif.crps) %>%
    group_by(names) %>%
    dplyr::mutate(mean = mean(dif.crps)) %>%
    ungroup()

  gg <- ggplot2::ggplot(t.crps.dif, ggplot2::aes(x = dif.crps, fill = names, color = names)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.4) +
    ggplot2::geom_vline(data = t.crps.dif, ggplot2::aes(xintercept = mean, color = names)) +
    ggplot2::theme(text = ggplot2::element_text(size = 20))
  ggplot2::ggsave(paste0(getwd(), filename, ".pdf"), gg, limitsize = FALSE)
}

#' This method filters the outcome of the [sim_e_values()] for the usual forecasts.
#' @param dt is the result of the [sim_e_values(usual.forecasts = TRUE)]. But it has to have the usual.forecasts = TRUE.
#' @export
filter_tibble_for_usual_forecasts <- function(dt) {
  return(dt %>%
           filter(!grepl("perfect-", names.F) & !grepl("perfect-", names.G)))
}

#' This method filters the outcome of the [sim_e_values()] for the usual forecasts and relocates the p-value-column to the end.
#' @param dt is the result of the [sim_e_values(usual.forecasts = TRUE)]. But it has to have the usual.forecasts = TRUE.
#' @export
filter_tibble_for_usual_forecasts_and_prettify <- function(dt) {
  return(filter_tibble_for_usual_forecasts(dt) %>%
           relocate(p.value.H0.rej, .after = last_col()))
}

#' This method prints the plots of a method and saves them to a file.
#' @param fileName of the resulting file.
#' @param fun = \(x) {function to be called (x) }, can be any function that produces an output.
#' @param data is the data to be put into fun. Usually return of [sim_e_values()].
#' @param paper = "A4", can be replaced with "A4r" for rotated A4.
#' @examples \dontrun{printPlot("test", \(x) {print_crps_diff_histogram_usual_forecasts(x)}, dt)}
#' @export
printPlot <- function(fileName, fun, data, paper = "A4") {
  pdf(paste0(getwd(), "/target/", fileName, ".pdf"), paper = paper)
  fun(data)
  dev.off()
}

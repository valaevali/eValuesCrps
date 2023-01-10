# Title     : Helper Scripts
# Created by: Velerie Haftka
# Created on: 26.01.2022

#' @export
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
  if (!("GRAPA" %in% method) &&
    !("lambda" %in% method) &&
    !("alternative" %in% method) &&
    !("alternative-mean" %in% method)) {
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

#' @export
check_input <- function(y, crps.F.para, crps.G.para, idx, method, lambda, p.value.method) {
  if (all(is.na(y)) || !is.vector(y)) {
    stop("The parameter y had to be a vector.")
  }

  check_input_crps_para(crps.F.para, "F")
  check_input_crps_para(crps.G.para, "G")

  if (is.na(idx)) {
    warning("Setting iteration to '1'")
    idx <- 1
  }
  if (!(any(c("lambda", "GRAPA", "alt-conf", "alt-cons", "alt-more-cons", "alt-mean") %in% method))) {
    warning("Parameter 'method' did not match any defined methods, setting to default 'lambda' for a fixed value")
    method <- list("lambda")
  }

  if (is.na(lambda) && ("lambda" %in% method)) {
    warning("Method with fixed lambda is choosen, but no lambda is provided. Setting to default lambda = 0.5")
    lambda <- 0.5
  }

  if (is.na(p.value.method)) {
    logger::log_debug("Parameter 'p.value.method' is na, hence no p.value will be calculated")
  }

  if (!is.na(p.value.method) &&
    !("t" %in% p.value.method) &&
    !("dm" %in% p.value.method)) {
    warning("Parameter 'p.value.method' has to be one of ('t','dm'). Setting to default 't'")
    p.value.method <- "t"
  }

  return(list("idx" = idx, "method" = method, "lambda" = lambda, "p.value.method" = p.value.method))
}

#' @export
check_input_crps_para <- function(crps.para, f.g) {
  if (all(is.na(crps.para))) {
    stop(sprintf("The CRPS parameter should not be empty for %s.", f.g))
  }
  if (!("mu" %in% names(crps.para)) &&
    !("points.cdf" %in% names(crps.para)) &&
    !("sd" %in% names(crps.para))) {
    stop(sprintf("The CRPS parameter ('mu','sd') or ('points.cdf') should not be empty for %s.", f.g))
  }
  if (("mu" %in% names(crps.para)) && ("sd" %in% names(crps.para))) {

    if (!is.numeric(crps.para$mu) && !is.vector(crps.para$mu)) {
      stop(sprintf("The CRPS parameter 'mu' should be a constant or a vector for %s.", f.g))
    }
    if (!is.numeric(crps.para$sd) && !is.vector(crps.para$sd)) {
      stop(sprintf("The CRPS parameter 'sd' should be a constant or a vector for %s.", f.g))
    }
    if ("w" %in% names(crps.para) &&
      !all(is.na(crps.para$w)) &&
      !is.numeric(crps.para$w) &&
      !(is.vector(crps.para$w) && !is.matrix(crps.para$w))) {
      stop(sprintf("The CRPS parameter 'w' should be either NA or a constant or a vector or a matrix for %s.", f.g))
    }

    if ((is.matrix(crps.para$mu) && (!is.matrix(crps.para$sd) || ("w" %in% names(crps.para) && !is.matrix(crps.para$w)))) ||
      (is.matrix(crps.para$sd) && (!is.matrix(crps.para$mu) || ("w" %in% names(crps.para) && !is.matrix(crps.para$w)))) ||
      (("w" %in% names(crps.para) && !is.matrix(crps.para$w)) && (!is.matrix(crps.para$sd) || !is.matrix(crps.para$mu)))) {
      stop(sprintf("One of the parameter for '%s' is a matrix, but not the others. Make sure all the parameters for one forecast have the same form", f.g))
    }
  } else {
    if (!is.numeric(crps.para$points.cdf) &&
      !is.vector(crps.para$points.cdf) &&
      all(is.na(crps.para$points.cdf$points)) &&
      all(is.na(crps.para$points.cdf$cdf))) {
      stop(sprintf("The CRPS parameter 'points.cdf' should be a matrix containing points and cdf for %s.", f.g))
    }
  }
}

#' Creates list for one input forecast of the form N(mu,sd)
#' @description This is a helper function which returns the correct form for the forecast input for [sim_e_values()].
#'
#' @param forecast.name is the name of the forecast.
#' @param mu is the mean of the forecast, can be constant, vector or matrix.
#' @param sd is the variance of the forecast, can be a constant, vector or matrix.
#' @param w = NA, optional, is the weight of a mixed forecast, can be a constant, vector or matrix.
#' @export
forecast_input <- function(mu, sd, w = NA) {
  if (all(is.na(w))) {
    return(list("mu" = mu, "sd" = sd))
  }
  return(list("mu" = mu, "sd" = sd, "w" = w))
}

#' Calculating Infimum of difference of two forecasts
#' @description This method makes a gitter search for global minima
#' @export
optim_inf_value <- function(f, start.points = 2, min.value = 0.0001, max.value = 1) {
  min(sapply(seq(min.value, max.value, length.out = start.points), \(y) { optim(y, f, lower = min.value, upper = max.value, method = "L-BFGS-B")$value }))
}

#' Creates helper functions for different forecasts
#' @description This methods creates the additional functions from the input parameters for the calculations of lambda.
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
create_crps_fun <- function(n.obs = 200, mu = 0, sd = 1, w = 1, points.cdf = NA, ...) {
  if (is.matrix(mu) || is.matrix(sd) || is.matrix(w)) {
    method <- 'mixnorm'
    crps.fun <- \(y) { scoringRules::crps_mixnorm(y = y, m = mu, s = sd, w = w) }
    crps.fun.y.matrix <- \(y) { sapply(1:dim(y)[2], \(i) { crps.fun(y[, i]) }) }
    sample.fun <- \(n) { matrix(rnorm(n * n.obs, mean = mu, sd = sd), nrow = n.obs) }
    inf.crps.fun <- \(x, j) { scoringRules::crps_mixnorm(y = x, m = as.matrix(t(mu[j,]), nrow = 1), s = as.matrix(t(sd[1,])), w = as.matrix(t(w[1,]))) }
  } else if (!all(is.na(points.cdf))) {
    method <- 'raw'
    crps.fun <- \(y) { sapply(seq_along(y), \(i) { crps_rf(y = y[i], points.cdf = points.cdf[i][[1]]) }) }

    crps.fun.y.matrix <- \(y) {
      result <- matrix(nrow = dim(y)[1], ncol = dim(y)[2])
      for (k in 1:dim(y)[1]) {
        for (l in 1:dim(y)[2]) {
          fun.y <- \(i, j) { crps_rf(y = y[i, j], points.cdf = points.cdf[j][[1]]) }
          result[k, l] <- fun.y(k, l)
        }
      }
      return(result)
    }

    sample.fun <- \(n, j) { rcdf_rf(points.cdf = points.cdf[j][[1]], n = n) }
    inf.crps.fun <- \(x, j) { crps_rf(y = x, points.cdf = points.cdf[j][[1]]) }
  } else {
    method <- 'norm'
    crps.fun <- \(y) { scoringRules::crps_norm(y = y, mean = mu, sd = sd) }
    crps.fun.y.matrix <- crps.fun
    sample.fun <- \(n) { matrix(rnorm(n * n.obs, mean = mu, sd = sd), nrow = n.obs) }

    if (length(mu) > 1 & length(sd) > 1) {
      inf.crps.fun <- \(x, j) { scoringRules::crps_norm(y = x, mean = mu[j], sd = sd[j]) }
    } else if (length(mu) > 1 & length(sd) == 1) {
      inf.crps.fun <- \(x, j) { scoringRules::crps_norm(y = x, mean = mu[j], sd = sd) }
    } else if (length(mu) == 1 & length(sd) > 1) {
      inf.crps.fun <- \(x, j) { scoringRules::crps_norm(y = x, mean = mu, sd = sd[j]) }
    } else {
      inf.crps.fun <- \(x, j) { scoringRules::crps_norm(y = x, mean = mu, sd = sd) }
    }
  }
  return(list("method" = method, "fun" = crps.fun, "crps.fun.y.matrix" = crps.fun.y.matrix,
              "sample.fun" = sample.fun, "inf.fun" = inf.crps.fun))
}

#' @export
crps_rf <- function(y, points.cdf) {
  # Check input
  if (!is.vector(y, "numeric"))
    stop("obs must be a numeric vector")
  if (length(y) != 1 && length(y) != length(points.cdf$points))
    stop("y must have length 1 or the same length as the predictions")

  w <- c(0, diff(points.cdf$cdf))
  a <- points.cdf$cdf + 0.5 * w
  crps0 <- function(y) 2 * sum(points.cdf$cdf * ((y < points.cdf$points) - a) * (points.cdf$points - y))
  sapply(y, crps0)
}

#' @export
rcdf_rf <- function(points.cdf, n) {
  r <- runif(n, min = min(points.cdf$points), max = max(points.cdf$points))

  cdf0 <- function(r) {
    # Evaluate CDF (stepfun) at thresholds
    stats::stepfun(x = points.cdf$points, y = c(0, points.cdf$cdf))(r)
  }

  sapply(r, cdf0)
}
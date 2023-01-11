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
#'
#' @description This is a helper function which returns the correct form for the forecast input, which are normaly
#' distributed, mixed normally distributed or are raw forecasts.
#'
#' @param mu is the mean of the forecast, can be constant, vector or matrix.
#' @param sd is the variance of the forecast, can be a constant, vector or matrix.
#' @param w = NA, optional, is the weight of a mixed forecast, can be a constant, vector or matrix.
#' @param points = NA, optional, are the points used for raw forecasts, see \code{\link{crps_rf}} and \code{\link{rcdf_rf}} for more details.
#' @param cdf = NA, optional, are the points used for raw forecasts, see \code{\link{crps_rf}} and \code{\link{rcdf_rf}} for more details.
#'
#' @return
#' Returns the list used in \code{\link{e_value}}.
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
forecast_input <- function(mu, sd, w = NA, points = NA, cdf = NA) {
  if (!all(is.na(points)) && !all(is.na(cdf))) {
    return(list("points.cdf" = tibble("points" = points, "cdf" = cdf)))
  }
  if (all(is.na(w))) {
    return(list("mu" = mu, "sd" = sd))
  }
  return(list("mu" = mu, "sd" = sd, "w" = w))
}

#' Calculating Infimum of difference of two forecasts
#'
#' @description This method makes a gitter search for global minima. The start points are taken in between the
#' \code{min.value} and the \code{max.value}.
#'
#' @param f the function which has to be optimized. The function has to take one value, ie f <- \(y) {}.
#' @param start.points = 2, how many points to start the grid search with
#' @param min.value = 0.0001, minimal value for the start.points to choose from
#' @param max.value = 1, maximum value for the start.points to choose from
#'
#' @return
#' Returns a single infimum value.
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
optim_inf_value <- function(f, start.points = 2, min.value = 0.0001, max.value = 1) {
  min(sapply(seq(min.value, max.value, length.out = start.points), \(y) { optim(y, f, lower = min.value, upper = max.value, method = "L-BFGS-B")$value }))
}

#' Creates helper functions for different forecasts
#'
#' @description This methods creates the additional functions from the input parameters for the calculations of lambda.
#'
#' @param n.obs = NA, number of observations, if not provided, this method tries to guess the n.obs from the other parameters provided.
#' @param mu = 0, the mean of the forecast.
#' @param sd = 1, the variance of the forecast.
#' @param w = 1, the weight of a mixed forecast.
#' @param points.cdf = NA, the \code{data.frame} for raw forecasts.
#' @param ... Additional parameters which are ignored.
#'
#' @returns
#' A list of method = ("norm"|"mixnorm"|"raw"), fun = \(y) {scoringRules::crps_norm(y=y, mean=mu, sd=sd} or \(y) { scoringRules::crps_mixnorm(y = y, m = mu, s = sd, w = w) } for mixed norm,
#' crps.fun.y.matrix is a function only special if method = "mixnorm" to calculate the alternative betting adaptive to forecasts otherwise it is equal to fun,
#' rnorm is the function to calculate randomly normally distributed variables, inf.fun is the function to calculate the infimum for this forecast.
#' For the raw method, see \code{\link{crps_rf}} or \code{\link{rcdf_rf}}.
#'
#' @examples
#' create_crps_fun(n.obs = 200, mu = rnorm(200), sd = 1)
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
create_crps_fun <- function(n.obs = NA, mu = 0, sd = 1, w = 1, points.cdf = NA, ...) {
  if (is.na(n.obs)) {
    if (is.vector(mu) && length(mu) > 1) n.obs <- length(mu)
    else if (is.matrix(mu) && nrow(mu) > 1) n.obs <- nrow(mu)
    else if (is.vector(sd) && length(sd) > 1) n.obs <- length(sd)
    else if (is.matrix(sd) && nrow(sd) > 1) n.obs <- nrow(sd)
    else if (is.vector(w) && length(w) > 1) n.obs <- length(w)
    else if (is.matrix(w) && nrow(w) > 1) n.obs <- nrow(w)
    else if (is.matrix(points.cdf) && nrow(points.cdf) > 1) n.obs <- nrow(points.cdf)
    else if (is.vector(mu) && length(mu) == 1 && is.vector(sd) && length(sd) == 1) n.obs <- 1
    else {
      stop("Number of rows cannot be determined! Please set n.obs.")
    }
  }

  if (is.matrix(mu) || is.matrix(sd) || is.matrix(w)) {
    method <- 'mixnorm'
    crps.fun <- \(y) { scoringRules::crps_mixnorm(y = y, m = mu, s = sd, w = w) }
    crps.fun.y.matrix <- \(y) { sapply(1:dim(y)[2], \(i) {scoringRules::crps_mixnorm(y = y[, i], m = mu, s = sd, w = w)}) }
    sample.fun <- \(n) { matrix(rnorm(n * n.obs, mean = mu, sd = sd), nrow = n.obs) }
    inf.crps.fun <- \(x, j) { scoringRules::crps_mixnorm(y = x, m = as.matrix(t(mu[j,]), nrow = 1), s = as.matrix(t(sd[1,])), w = as.matrix(t(w[1,]))) }
  } else if (!all(is.na(points.cdf))) {
    method <- 'raw'
    crps.fun <- \(y) { sapply(seq_along(y), \(i) { crps_rf(y = y[i], points.cdf = points.cdf[i][[1]]) }) }

    crps.fun.y.matrix <- \(y) {
      if (all(is.na(dim(y)))) dim <- c(1, length(y)) else dim <- dim(y)
      result <- matrix(nrow = dim[1], ncol = dim[2])
      for (k in 1:dim[1]) {
        for (l in 1:dim[2]) {
          if (all(is.na(dim(y)))) {
            fun.y <- \(i, j) { crps_rf(y = y[j], points.cdf = points.cdf[i][[1]]) }
          } else {
            fun.y <- \(i, j) { crps_rf(y = y[i, j], points.cdf = points.cdf[i][[1]]) }
          }
          result[k,l] <- fun.y(k, l)
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

#' Continuous ranked probability score (CRPS)
#'
#' @description Computs the CRPS of raw forecasts.
#'
#'
#' @param y a numeric vector of observations of the same length as the number of points.cdf, or of length 1.
#' @param points.cdf a \code{data.frame} of numeric variables, used to compute the empirical distribution of the
#' variables in \code{points.cdf}.
#'
#' @details
#' This function uses adapted code taken from the function \code{crps_edf} of the \pkg{scoringRules} package and of the
#' function \code{cdf} of the \pkg{isoditrreg} package.
#'
#' @return
#' If the input y is only a scalar, then it does return a single CRPS value. If y is a vector, it does return a vector
#' of CRPS values, evaluated for each y.
#'
#' @seealso \code{\link{e_value}}
#'
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

#' Random values of the cumulative distribution function (CDF) of raw forecasts
#'
#' @description Evaluate the cumulative distribution function (CDF) of raw forecasts in a \code{data.frame} at
#' randomly generated thresholds.
#'
#'
#' @param points.cdf a \code{data.frame} of numeric variables, used to compute the empirical distribution of the
#' variables in \code{points.cdf}.
#' @param n is the number of randomly generated thresholds to evaluate the cdf at.
#'
#' @details
#' The CDFs are considered as piecewise constant stepfunctions. The \code{points} in the \code{data.frame}
#' \code{points.cdf} are the points where the empirical distribution of the forecasts has jumps and \code{cdf} in the
#' \code{data.frame} \code{points.cdf} are the corresponding CDF values.
#'
#' @return
#' A vector of probabilities giving the evaluated CDFs at the randomly generated thresholds.
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
rcdf_rf <- function(points.cdf, n) {
  r <- runif(n, min = min(points.cdf$points), max = max(points.cdf$points))

  cdf0 <- function(r) {
    # Evaluate CDF (stepfun) at thresholds
    stats::stepfun(x = points.cdf$points, y = c(0, points.cdf$cdf))(r)
  }

  sapply(r, cdf0)
}
#' E-value
#'
#' @description
#' This method is the main method of this package and calculates the e.values for one input.
#' This method uses the optional stopping principle, because of \code{max(cumprod(e.value)}. If your goal is to calculate
#' multiple e-values sequentially, then run this method once and for the next steps, use \code{\link{next_k_e_values_for_point_cdfs}},
#' because it does optimize the computational costs massively.
#'
#' @usage e_value(y, crps.F.para, crps.G.para)
#'
#' @param y oberservations y \in R^k
#' @param crps.F.para of the form \code{list("mu" = mu, "sd" = 1)} or of the form \code{list("points.cdf" = tibble with points and cdfs)} (see \code{\link{crps_rf}} and \code{\link{rcdf_rf}} for more details on raw forecasts)
#' @param crps.G.para of the form \code{list("mu" = mu, "sd" = 1)} or of the form \code{list("points.cdf" = tibble with points and cdfs)} (see \code{\link{crps_rf}} and \code{\link{rcdf_rf}} for more details on raw forecasts)
#' @param idx = 1, this is an aditional parameter to idex the run of this call
#' @param method = list("GRAPA", "lambda", "alt-conf", "alt-cons", "alt-more-cons"), is a list containing all the method names for calculating the different lambdas
#' @param lambda = 0.5, lambda entry for a fixed value
#' @param p.value.method = NA, can be "t" for t-test and "dm" for dm.test of the package forecast
#' @param old.run.e.value is the return of the last call to \code{e_value}, this can be called directly or with the function
#' \code{\link{next_k_e_values_for_point_cdfs}}.
#' @param k is the parameter used for \code{\link{next_k_e_values_for_point_cdfs}} to specify how many new observations are evaulated
#'
#' @examples
#' mu <- rnorm(10)
#' y <- rnorm(10, mu)
#' e_value(y = y, crps.F.para = list("mu" = mu, "sd" = 1), crps.G.para = list("mu" = 0, "sd" = 2))
#'
#' @return
#' Returns a list containing the input values and the calculated e-values and p-values (if specified).
#'
#' @seealso \code{\link{next_k_e_values_for_point_cdfs}}, \code{\link{crps_rf}}, \code{\link{rcdf_rf}}
#'
#' @export
e_value <- function(y, crps.F.para, crps.G.para, idx = 1,
                    method = list("alt-cons"), lambda = 0.5, p.value.method = NA, old.run.e.value = NA, k = NA) {
  checkedInput <- check_input(y, crps.F.para, crps.G.para, idx, method, lambda, p.value.method)
  idx <- checkedInput$idx
  method <- checkedInput$method
  lambda <- checkedInput$lambda
  p.value.method <- checkedInput$p.value.method

  crps.F.para <- base::append(crps.F.para, rlang::exec(create_crps_fun, length(y), !!!crps.F.para))
  crps.G.para <- base::append(crps.G.para, rlang::exec(create_crps_fun, length(y), !!!crps.G.para))
  n.obs <- length(y)

  if (!is.na(k)) {
    crps.F <- c(old.run.e.value$crps.F, crps.F.para$fun(y[(n.obs - k + 1):n.obs]))
    crps.G <- c(old.run.e.value$crps.G, crps.G.para$fun(y[(n.obs - k + 1):n.obs]))
  } else {
    crps.F <- crps.F.para$fun(y)
    crps.G <- crps.G.para$fun(y)
  }

  # Calculating inf.crps
  logger::log_debug("Starting infimum caclulation")
  inf.crps <- get_inf_crps(crps.F.para, crps.G.para, n.obs, k, if (!all(is.na(old.run.e.value))) old.run.e.value$inf.crps else NA)

  T.F.G <- (crps.F - crps.G) / inf.crps
  e.values <- list("crps.F.fun" = crps.F.para, "crps.F" = crps.F, "crps.G.fun" = crps.G.para, "crps.G" = crps.G,
                   "inf.crps" = inf.crps, "y" = y, "idx" = idx, "lambda" = lambda, "method" = method, "p.value.method" = p.value.method)

  logger::log_debug("Starting lambda caclulation")
  if ("lambda" %in% method) {
    # Lambda is fix
    e.value <- 1 + lambda * T.F.G
    e.value.prod <- max(cumprod(e.value))
    e.values <- base::append(e.values, list("e.value.lambda" = e.value, "e.value.lambda.prod" = e.value.prod))
  }

  logger::log_debug("Starting GRAPA caclulation")
  if ("GRAPA" %in% method) {
    # GRAPA
    e.values <- base::append(e.values, e_value_calculate_lambda_for_grapa_betting(T.F.G))
  }

  logger::log_debug("Starting alternative betting caclulation")
  if (any(c("alt-conf", "alt-cons", "alt-more-cons") %in% method)) {
    # Alternative
    e.values <- base::append(e.values, e_value_calculate_lambda_for_alternative_betting(T.F.G, crps.F.para, crps.G.para, inf.crps, method, old.run.e.value, k))
  }

  logger::log_debug("Starting p-value caclulation")
  if (!is.na(p.value.method)) {
    # P-value
    e.values <- base::append(e.values, p_value_t_test(crps.F, crps.G, p.value.method))
  }

  return(e.values)
}

#' E-value
#'
#' @description
#' This method prepares the information of the old run from \code{\link{e_value}} to be run again and runs it with some additional input.
#' With this approach, all the infimums, do not have to be calculated for each k, but only for the new ones and also for the betting
#' approach, not the whole matrices have to be calculated again. This procedure reduces the computational cost massively.
#'
#' @usage next_k_e_values_for_point_cdfs(e.value.run.before, new.y, new.crps.F.para, new.crps.G.para)
#'
#' @param e.value.run.before is the result of the run of \code{\link{e_value}}.
#' @param new.y oberservations y \in R^k, is the new observation to evaluate
#' @param new.crps.F.para of the form \code{list("mu" = mu[k], "sd" = NA)} or of the form
#' \code{list("points.cdf" = tibble with points and cdfs)}, if either "mu" or "sd" are a single value, set it to NA.
#' @param new.crps.G.para of the form \code{list("mu" = mu[k], "sd" = NA)} or of the form
#' \code{list("points.cdf" = tibble with points and cdfs)}, if either "mu" or "sd" are a single value, set it to NA.
#' @param idx = 2, this is an aditional parameter to idex the run of this call
#'
#' @examples
#' mu <- rnorm(10)
#' y <- rnorm(10, mu)
#' result <- e_value(y = y, crps.F.para = list("mu" = mu, "sd" = 1), crps.G.para = list("mu" = 0, "sd" = 2))
#' new.mu <- rnorm(1)
#' next_k_e_values(result, new.y = rnorm(1, new.mu), cprs.F.para = list("mu" = new.mu, "sd" = NA), crps.G.para = list("mu" = NA, "sd" = NA))
#'
#' @return
#' Returns a list containing the input values and the calculated e-values.
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
next_k_e_values_for_point_cdfs <- function(e.value.run.before, new.y, new.crps.F.para, new.crps.G.para, idx = 2) {
  if (all(is.na(e.value.run.before))) {
    stop("Initial step needs the return of the function e_values as input!")
  }
  if (all(is.na(new.y))) {
    stop("new.y cannot be NA, please provide new.y=?")
  }

  k <- length(new.y)

  y <- e.value.run.before$y
  crps.F.para <- e.value.run.before$crps.F.fun
  crps.G.para <- e.value.run.before$crps.G.fun
  method <- e.value.run.before$method[[1]]
  lambda <- e.value.run.before$lambda
  p.value.method <- e.value.run.before$p.value.method

  y <- append(y, new.y)
  if (("norm" == crps.F.para$method) && ("mu" %in% names(crps.F.para)) && ("sd" %in% names(crps.F.para))) {
    crps.F.para <- list("mu" = na.omit(append(crps.F.para$mu, new.crps.F.para$mu)), "sd" = na.omit(append(crps.F.para$mu, new.crps.F.para$mu)))
  } else if ("points.cdf" %in% names(crps.F.para) && "raw" == crps.F.para$method) {
    crps.F.para <- list("points.cdf" = append(crps.F.para$points.cdf, new.crps.F.para))
  } else {

  }
  if (("norm" == crps.G.para$method) && ("mu" %in% names(crps.G.para)) && ("sd" %in% names(crps.G.para))) {
    crps.G.para <- list("mu" = na.omit(append(crps.G.para$mu, new.crps.G.para$mu)), "sd" = na.omit(append(crps.G.para$mu, new.crps.G.para$mu)))
  } else if ("points.cdf" %in% names(crps.G.para) && "raw" == crps.G.para$method) {
    crps.G.para <- list("points.cdf" = append(crps.G.para$points.cdf, new.crps.G.para))
  } else {

  }

  return(e_value(y = y, crps.F.para = crps.F.para, crps.G.para = crps.G.para, idx = idx,
                 method = method, lambda = lambda, p.value.method = p.value.method, old.run.e.value = e.value.run.before, k = k))
}

#' Infimum calculation
#'
#' @description
#' This method calculates the infimum of crps.F - crps.G.
#'
#' @param crps.F.para of the form \code{list("mu" = mu, "sd" = 1)} or \code{list("points.cdf" = tibble with points and cdfs)}
#' @param crps.G.para of the form \code{ist("mu" = -mu, "sd" = 1)} or \code{list("points.cdf" = tibble with points and cdfs)}
#' @param n.obs is the number of observations
#' @param k is the number of new observations for \code{\link{next_k_e_values_for_point_cdfs}}
#' @param old.inf is the old infimum value of the last run of \code{e_value()}, used by \code{\link{next_k_e_values_for_point_cdfs}}
#'
#' @return
#' Returns the infimum value for crps.F - crps.G.
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
get_inf_crps <- function(crps.F.para, crps.G.para, n.obs, k = NA, old.inf = NA) {
  # Check if F and G are of the form N(mu,sigma^2)
  if (crps.F.para$method == 'norm' &&
    crps.G.para$method == 'norm') {
    # Check if F.sd != G.sd
    if (crps.F.para$sd != crps.G.para$sd) {
      # Here we have to take the minimum of the null.value vector => otherwise this is not the minimum
      null.value <- min((sqrt(2 * crps.F.para$sd) * crps.G.para$mu - sqrt(2 * crps.G.para$sd) * crps.F.para$mu) /
                          (sqrt(2 * crps.F.para$sd) - sqrt(2 * crps.G.para$sd)))
      return(abs(min(crps.F.para$fun(null.value) - crps.G.para$fun(null.value))))
    } else {
      return(abs(min(crps.F.para$mu - crps.G.para$mu)))
    }
  } else {

    optim.inf.fun <- \(i) { optim_inf_value(\(x) { crps.F.para$inf.fun(x, i) - crps.G.para$inf.fun(x, i) },
                                            min.value = -10, max.value = 10) }

    if (!is.na(k) & !is.na(old.inf)) {
      new.inf <- for(i in (n.obs - k + 1):n.obs) optim.inf.fun(i)
      return(abs(min(c(old.inf, new.inf))))
    }
    return(abs(min(sapply((1:n.obs), optim.inf.fun))))
  }
}

#' Betting strategies
#'
#' @description
#' This method calculates the lambda with the grapa betting strategy
#'
#' @param T.F.G = (crps.F - crps.G) / inf.crps
#'
#' @return
#' Returns a list containing the calculated \code{e.value.grapa}, \code{e.value.grapa.prod} and the corresponding
#' \code{lambda.grapa}. \code{e.value.grapa} and \code{lambda.grapa} are both vectors of the size of the observation.
#' \code{e.value.grapa.prod} is the e-value process at time k (number of observations).
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
e_value_calculate_lambda_for_grapa_betting <- function(T.F.G) {
  t.F.G.i <- cumsum(T.F.G)
  t.F.G.i.2 <- cumsum((T.F.G)^2)

  lambda.grapa <- t.F.G.i / (t.F.G.i.2)

  # Here we have to manually add for lambda[1:10] <- 0.5, since the sum should only take [t-1] and flatten it
  lambda.grapa <- append(0.5, lambda.grapa)[-(length(lambda.grapa) + 1)]
  lambda.grapa[2:10] <- 0.5
  lambda.grapa[which(lambda.grapa < 0.0001)] <- 0.0001
  lambda.grapa[which(lambda.grapa > 1)] <- 1

  e.value.grapa <- 1 + lambda.grapa * T.F.G
  e.value.prod.grapa <- max(cumprod(e.value.grapa))

  return(list("e.value.grapa" = e.value.grapa, "e.value.grapa.prod" = e.value.prod.grapa, "lambda.grapa" = lambda.grapa))
}

#' Betting strategies
#'
#' @description
#' This method calculates the lambda with the alternative betting strategy
#'
#' @param T.F.G = (crps.F - crps.G) / inf.crps
#' @param crps.F.para of the form \code{list("mu" = mu, "sd" = 1)} or \code{list("points.cdf" = tibble with points and cdfs)}
#' @param crps.G.para of the form \code{ist("mu" = -mu, "sd" = 1)} or \code{list("points.cdf" = tibble with points and cdfs)}
#' @param inf.crps = infimum of crps.F - crps.G over y
#' @param method = list("alt-conf", "alt-cons", "alt-more-cons"), is a list containing all the method names for calculating the
#' different lambdas
#' @param old.run.e.value is the return of the last call to \code{\link{e_value}}, this can be called directly or with the function
#' \code{\link{next_k_e_values_for_point_cdfs}}
#' @param k is the parameter used for \code{\link{next_k_e_values_for_point_cdfs}} to specify how many new observations are evaulated
#'
#' @return
#' Returns a list for each alternative betting strategy specified, containing \code{e.value.alt.(conf|cons|more.cons)},
#' \code{e.value.alt.(conf|cons|more.cons).prod} and \code{lambda.alt.(conf|cons|more.cons)}.
#' \code{e.value.alt.(conf|cons|more.cons)} and \code{lambda.alt.(conf|cons|more.cons)} are both vectors of the size of
#' the observation. \code{e.value.alt.(conf|cons|more.cons).prod} is the e-value process at time k (number of observations).
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
e_value_calculate_lambda_for_alternative_betting <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, method, old.run.e.value, k) {
  result <- NA

  if ("alt-conf" %in% method) {
    result <-
      base::append(result,
                   e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps,
                                                                         suffix = "alt.conf", F.proportion = 0,
                                                                         G.proportion = 1,
                                                                         crps.alt.old = if (!all(is.na(old.run.e.value))) old.run.e.value$crps.alt.conf else NA, k = k))
  }

  if ("alt-cons" %in% method) {
    result <-
      base::append(result,
                   e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps,
                                                                         suffix = "alt.cons", F.proportion = 0.15,
                                                                         G.proportion = 0.85,
                                                                         crps.alt.old = if (!all(is.na(old.run.e.value))) old.run.e.value$crps.alt.cons else NA, k = k))
  }

  if ("alt-more-cons" %in% method) {
    result <-
      base::append(result,
                   e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps,
                                                                         suffix = "alt.more.cons", F.proportion = 0.25,
                                                                         G.proportion = 0.75,
                                                                         crps.alt.old = if (!all(is.na(old.run.e.value))) old.run.e.value$crps.alt.more.cons else NA, k = k))
  }

  result <- result[!is.na(result)]
  return(result)
}

e_value_calculate_lambda_for_alternative_betting_each <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, suffix, F.proportion, G.proportion, crps.alt.old, k) {
  min.sample <- 50

  if (!is.na(k) & !all(is.na(crps.alt.old)) && crps.F.para$method == 'raw') {
    y.sim <- (F.proportion * t(sapply((length(crps.F.para$points.cdf) - k + 1):length(crps.F.para$points.cdf), \(i) { crps.F.para$sample.fun(min.sample, i) })) +
      G.proportion * t(sapply((length(crps.F.para$points.cdf) - k + 1):length(crps.G.para$points.cdf), \(i) { crps.G.para$sample.fun(min.sample, i) })))
    crps.alt.new <- crps.F.para$crps.fun.y.matrix(y.sim) - crps.G.para$crps.fun.y.matrix(y.sim)
    crps.alt <- rbind(crps.alt.old, crps.alt.new)
  } else {
    if (crps.F.para$method == 'raw') {
      y.sim.F <- t(sapply(seq_along(crps.F.para$points.cdf), \(i) { crps.F.para$sample.fun(min.sample, i) }))
    } else {
      y.sim.F <- crps.F.para$sample.fun(min.sample)
    }
    if (crps.G.para$method == 'raw') {
      y.sim.G <- t(sapply(seq_along(crps.G.para$points.cdf), \(i) { crps.G.para$sample.fun(min.sample, i) }))
    } else {
      y.sim.G <- crps.G.para$sample.fun(min.sample)
    }

    y.sim <- (F.proportion * y.sim.F + G.proportion * y.sim.G)
    crps.alt <- crps.F.para$crps.fun.y.matrix(y.sim) - crps.G.para$crps.fun.y.matrix(y.sim)
  }

  if (is.vector(crps.alt)) {
    lambda.alt <- mean(crps.alt / inf.crps) / mean((crps.alt / inf.crps)^2)
  } else {
    lambda.alt <- rowMeans(crps.alt / inf.crps) / rowMeans((crps.alt / inf.crps)^2)
  }
  lambda.alt[which(lambda.alt < 0.0001)] <- 0.0001
  lambda.alt[which(lambda.alt > 1)] <- 1

  e.value.alt <- 1 + lambda.alt * T.F.G
  e.value.alt.prod <- max(cumprod(e.value.alt))
  result <- setNames(list(e.value.alt, e.value.alt.prod, lambda.alt, crps.alt), c(paste0("e.value.", suffix), paste0("e.value.", suffix, ".prod"), paste0("lambda.", suffix), paste0("crps.", suffix)))
  return(result)
}

#' P-value calculations
#'
#' @description
#' This method calculates the Diebold-Mariano t-test
#'
#' @param crps.F the crps vector of the F forecast
#' @param crps.G the crps vector of the G forecast
#' @param p.value.method = one of ("t","dm"), "t" stand for t-test directly and dm stands for the dm.test method of the forecast package
#'
#' @return
#' Returns a single p-value.
#'
#' @seealso \code{\link{e_value}}
#'
#' @export
p_value_t_test <- function(crps.F, crps.G, p.value.method = "t") {
  if (length(crps.F) == 1 || length(crps.G) == 1) {
    stop("P-value cannot be calculated for one observation, please provide a larger data set.")
  }
  if (p.value.method == "dm" || !(is.numeric(crps.F) && is.numeric(crps.G))) {
    p.value <- as.numeric(forecast::dm.test(crps.F, crps.G, alternative = "greater")$p.value)
  } else {
    sigma.n <- 1 / length(crps.F) * sum((crps.F - crps.G)^2)
    test.statistic <- sqrt(length(crps.F)) * ((crps.F - crps.G) / sigma.n)
    p.value <- as.numeric(t.test(test.statistic, alternative = "greater")$p.value)
  }
  return(list("p.value" = p.value))
}
#' This method is the main method of this package and calculates the e.values for one iterationstep.
#' This method uses the optional stopping principle, because of max(cumprod(e.value).
#'
#' @param y oberservations y \in R^k
#' @param mu is the mean parameter for y ~ N(mu, 1)
#' @param crps.F.para of the form list("mu" = mu, "sd" = 1), here mu and sd do not have to be of the same form, but for mixnorm, they all have to have the same dimensions
#' @param crps.G.para of the form list("mu" = mu, "sd" = 1), here mu and sd do not have to be of the same form, but for mixnorm, they all have to have the same dimensions
#' @param idx = 1, this is an aditional parameter to idex the run of this call
#' @param method = list("GRAPA", "lambda", "alt-conf", "alt-cons", "alt-more-cons"), is a list containing all the method names for calculating the different lambdas
#' @param lambda = 0.5, lambda entry for a fixed value
#' @param p.value.method = NA, can be "t" for t-test and "dm" for dm.test of the package forecast
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

#' This method calculates the infimum of crps.F - crps.G.
#' @param crps.F.para of the form list("mu" = mu, "sd" = 1)
#' @param crps.G.para of the form list("mu" = -mu, "sd" = 1)
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

#' This method calculates the lambda with the grapa betting strategy
#' @param T.F.G = (crps.F - crps.G) / inf.crps
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

#' This method calculates the lambda with the alternative betting strategy
#' @param T.F.G = (crps.F - crps.G) / inf.crps
#' @param crps.F.para of the form list("mu" = mu, "sd" = 1)
#' @param crps.G.para of the form list("mu" = -mu, "sd" = 1)
#' @param inf.crps = infimum of crps.F - crps.G over y
#' @param method = list("alt-conf", "alt-cons", "alt-more-cons"), is a list containing all the method names for calculating the different lambdas
#' @export
e_value_calculate_lambda_for_alternative_betting <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, method, input, k) {
  result <- NA

  if ("alt-conf" %in% method) {
    result <-
      base::append(result,
                   e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps,
                                                                         suffix = "alt.conf", F.proportion = 0,
                                                                         G.proportion = 1,
                                                                         crps.alt.old = if (!all(is.na(input))) input$crps.alt.conf else NA, k = k))
  }

  if ("alt-cons" %in% method) {
    result <-
      base::append(result,
                   e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps,
                                                                         suffix = "alt.cons", F.proportion = 0.15,
                                                                         G.proportion = 0.85,
                                                                         crps.alt.old = if (!all(is.na(input))) input$crps.alt.cons else NA, k = k))
  }

  if ("alt-more-cons" %in% method) {
    result <-
      base::append(result,
                   e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps,
                                                                         suffix = "alt.more.cons", F.proportion = 0.25,
                                                                         G.proportion = 0.75,
                                                                         crps.alt.old = if (!all(is.na(input))) input$crps.alt.more.cons else NA, k = k))
  }

  result <- result[!is.na(result)]
  return(result)
}

e_value_calculate_lambda_for_alternative_betting_each <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, suffix, F.proportion, G.proportion, crps.alt.old, k) {
  min.sample <- 50

  if (!is.na(k) & !all(is.na(crps.alt.old))) {
    y.sim <- (F.proportion * sapply((length(crps.F.para$points.cdf) - k + 1):length(crps.F.para$points.cdf), \(i) { crps.F.para$sample.fun(min.sample, i) })) +
      G.proportion * sapply((length(crps.F.para$points.cdf) - k + 1):length(crps.G.para$points.cdf), \(i) { crps.G.para$sample.fun(min.sample, i) })
    crps.alt.new <- crps.F.para$crps.fun.y.matrix(y.sim) - crps.G.para$crps.fun.y.matrix(y.sim)
    crps.alt <- rbind(crps.alt.old, t(crps.alt.new))
  } else {
    if (crps.F.para$method == 'raw') {
      y.sim.F <- sapply(seq_along(crps.F.para$points.cdf), \(i) { crps.F.para$sample.fun(min.sample, i) })
    } else {
      y.sim.F <- crps.F.para$sample.fun(min.sample)
    }
    if (crps.G.para$method == 'raw') {
      y.sim.G <- sapply(seq_along(crps.G.para$points.cdf), \(i) { crps.G.para$sample.fun(min.sample, i) })
    } else {
      y.sim.G <- crps.G.para$sample.fun(min.sample)
    }

    y.sim <- t((F.proportion * y.sim.F + G.proportion * y.sim.G))
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

#' This method calculates the Diebold-Mariano t-test
#' @param crps.F the crps vector of the F forecast
#' @param crps.G the crps vector of the G forecast
#' @param p.value.method = one of ("t","dm"), "t" stand for t-test directly and dm stands for the dm.test method of the forecast package
#' @export
p_value_t_test <- function(crps.F, crps.G, p.value.method = "t") {
  if (p.value.method == "dm" || !(is.numeric(crps.F) && is.numeric(crps.G))) {
    p.value <- as.numeric(forecast::dm.test(crps.F, crps.G, alternative = "greater")$p.value)
  } else {
    sigma.n <- 1 / length(crps.F) * sum((crps.F - crps.G)^2)
    test.statistic <- sqrt(length(crps.F)) * ((crps.F - crps.G) / sigma.n)
    p.value <- as.numeric(t.test(test.statistic, alternative = "greater")$p.value)
  }
  return(list("p.value" = p.value))
}

#' @export
next_k_e_values_for_point_cdfs <- function(e.value.run.before, new.y, new.crps.F.para, new.crps.G.para, idx = 2, k = 1) {
  if (all(is.na(e.value.run.before))) {
    stop("Initial step needs the return of the function e_values as input!")
  }
  if (all(is.na(k)) || k <= 0) {
    stop("K must be a positive integer!")
  }
  if (all(is.na(new.y))) {
    stop("next.y cannot be NA, please provide next.obs=?")
  }

  y <- e.value.run.before$y
  crps.F.para <- e.value.run.before$crps.F.fun
  crps.G.para <- e.value.run.before$crps.G.fun
  method <- e.value.run.before$method[[1]]
  lambda <- e.value.run.before$lambda
  p.value.method <- e.value.run.before$p.value.method

  y <- append(y, new.y)
  crps.F.para <- list("points.cdf" = append(crps.F.para$points.cdf, new.crps.F.para))
  crps.G.para <- list("points.cdf" = append(crps.G.para$points.cdf, new.crps.G.para))

  return(e_value(y = y, crps.F.para = crps.F.para, crps.G.para = crps.G.para, idx = idx,
                 method = method, lambda = lambda, p.value.method = p.value.method, old.run.e.value = e.value.run.before, k = k))
}
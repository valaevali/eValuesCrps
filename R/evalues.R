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
                    method = list("alt-cons"), lambda = 0.5, p.value.method = NA) {
  checkedInput <- check_input(y, crps.F.para, crps.G.para, idx, method, lambda, p.value.method)
  idx <- checkedInput$idx
  method <- checkedInput$method
  lambda <- checkedInput$lambda
  p.value.method <- checkedInput$p.value.method

  crps.F.para <- base::append(crps.F.para, rlang::exec(create_crps_fun, length(y), !!!crps.F.para))
  crps.G.para <- base::append(crps.G.para, rlang::exec(create_crps_fun, length(y), !!!crps.G.para))
  n.obs <- length(y)

  crps.F <- crps.F.para$fun(y)
  crps.G <- crps.G.para$fun(y)

  # Calculating inf.crps
  logger::log_debug("Starting infimum caclulation")
  inf.crps <- get_inf_crps(crps.F.para, crps.G.para, n.obs)

  T.F.G <- (crps.F - crps.G) / inf.crps
  e.values <- list("crps.F.fun" = crps.F.para, "crps.F" = crps.F, "crps.G.fun" = crps.G.para, "crps.G" = crps.G,
                   "inf.crps" = inf.crps, "y" = y, "idx" = idx, "lambda" = lambda)

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
    e.values <- base::append(e.values, e_value_calculate_lambda_for_alternative_betting(T.F.G, crps.F.para, crps.G.para, inf.crps, method))
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
get_inf_crps <- function(crps.F.para, crps.G.para, n.obs) {
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
    return(abs(min(sapply((1:n.obs), \(i) { optim_inf_value(\(x) { crps.F.para$inf.fun(x, i) - crps.G.para$inf.fun(x, i) },
                                                            min.value = -10, max.value = 10) }))))
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
e_value_calculate_lambda_for_alternative_betting <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, method) {
  result <- NA

  if ("alt-conf" %in% method) {
    result <- base::append(result, e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps, suffix = "alt.conf", F.proportion = 0, G.proportion = 1))
  }

  if ("alt-cons" %in% method) {
    result <- base::append(result, e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps, suffix = "alt.cons", F.proportion = 0.15, G.proportion = 0.85))
  }

  if ("alt-more-cons") {
    result <- base::append(result, e_value_calculate_lambda_for_alternative_betting_each(T.F.G = T.F.G, crps.F.para = crps.F.para,
                                                                                         crps.G.para = crps.G.para, inf.crps = inf.crps, suffix = "alt.more.cons", F.proportion = 0.25, G.proportion = 0.75))
  }

  result <- result[!is.na(result)]
  return(result)
}

e_value_calculate_lambda_for_alternative_betting_each <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, suffix, F.proportion, G.proportion) {
  min.sample <- if (length(T.F.G) == 1) 20 else length(T.F.G)

  y.sim <- (F.proportion * crps.F.para$sample.fun(min.sample) + G.proportion * crps.G.para$sample.fun(min.sample))
  crps.alt <- crps.F.para$crps.fun.y.matrix(y.sim) - crps.G.para$crps.fun.y.matrix(y.sim)
  if (is.vector(crps.alt)) {
    lambda.alt <- mean(crps.alt / inf.crps) / mean((crps.alt / inf.crps)^2)
  } else {
    lambda.alt <- rowMeans(crps.alt / inf.crps) / rowMeans((crps.alt / inf.crps)^2)
  }
  lambda.alt[which(lambda.alt < 0.0001)] <- 0.0001
  lambda.alt[which(lambda.alt > 1)] <- 1

  e.value.alt <- 1 + lambda.alt * T.F.G
  e.value.alt.prod <- max(cumprod(e.value.alt))
  return(setNames(list(e.value.alt, e.value.alt.prod, lambda.alt), c(paste0("e.value.", suffix), paste0("e.value.", suffix, ".prod"), paste0("lambda.", suffix))))
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

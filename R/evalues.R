#' This method is the main method of this package and calculates the e.values for one iterationstep.
#' This method uses the optional stopping principle, because of max(cumprod(e.value).
#'
#' @param y oberservations y \in R^k
#' @param mu is the mean parameter for y ~ N(mu, 1)
#' @param crps.F.para of the form list("mu" = mu, "sd" = 1) plus the list created with [create_crps_fun()]
#' @param crps.G.para of the form list("mu" = mu, "sd" = 1) plus the list created with [create_crps_fun()]
#' @param it = 1, is only importend when used with [sim_e_values()] and specifies which iteration it is
#' @param method = list("GRAPA", "lambda", "alternative", "alternative-mean"), is a list containing all the method names for calculating the different lambdas
#' @param lambda = 0.5, lambda entry for a fixed value
#' @param p.value.method = NA, can be "t" for t-test and "dm" for dm.test of the package forecast
#' @export
e_value <- function(y, crps.F.para, crps.G.para, it = 1,
                    method = list("GRAPA", "lambda", "alternative", "alternative-mean"), lambda = 0.5, p.value.method = NA) {
  checkedInput <- check_input(y, crps.F.para, crps.G.para, it, method, lambda, p.value.method)
  it <- checkedInput$it
  method <- checkedInput$method
  lambda <- checkedInput$lambda
  p.value.method <- checkedInput$p.value.method

  crps.F.para <- base::append(crps.F.para, rlang::exec(create_crps_fun, length(y), !!!crps.F.para))
  crps.G.para <- base::append(crps.G.para, rlang::exec(create_crps_fun, length(y), !!!crps.G.para))

  crps.F <- crps.F.para$fun(y)
  crps.G <- crps.G.para$fun(y)

  # Calculating inf.crps
  inf.crps <- get_inf_crps(crps.F.para, crps.G.para)

  T.F.G <- (crps.F - crps.G) / inf.crps
  e.values <- list("crps.F.fun" = crps.F.para, "crps.F" = crps.F, "crps.G.fun" = crps.G.para, "crps.G" = crps.G,
                   "inf.crps" = inf.crps, "y" = y, "it" = it, "lambda" = lambda)

  if ("lambda" %in% method) {
    # Lambda is fix
    e.value <- 1 + lambda * T.F.G
    e.value.prod <- 1 / max(cumprod(e.value))
    e.values <- base::append(e.values, list("e.value.lambda" = e.value, "e.value.lambda.prod" = e.value.prod))
  }

  if ("GRAPA" %in% method) {
    # GRAPA
    e.values <- base::append(e.values, e_value_calculate_lambda_for_grapa_betting(T.F.G))
  }

  if ("alternative" %in% method || "alternative-mean" %in% method) {
    # Alternative
    e.values <- base::append(e.values, e_value_calculate_lambda_for_alternative_betting(T.F.G, crps.F.para, crps.G.para, inf.crps, method))
  }

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
get_inf_crps <- function(crps.F.para, crps.G.para) {
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
    return(abs(optim_inf_value(\(x) { crps.F.para$inf.fun(x, 1) - crps.G.para$inf.fun(x, 1) },
                               min.value = -10, max.value = 10)))
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
  e.value.prod.grapa <- 1 / max(cumprod(e.value.grapa))

  return(list("e.value.grapa" = e.value.grapa, "e.value.grapa.prod" = e.value.prod.grapa, "lambda.grapa" = lambda.grapa))
}

#' This method calculates the lambda with the alternative betting strategy
#' @param T.F.G = (crps.F - crps.G) / inf.crps
#' @param crps.F.para of the form list("mu" = mu, "sd" = 1)
#' @param crps.G.para of the form list("mu" = -mu, "sd" = 1)
#' @param inf.crps = infimum of crps.F - crps.G over y
#' @param method = list("GRAPA", "lambda", "alternative", "alternative-mean"), is a list containing all the method names for calculating the different lambdas, only parameter used in this method is 'alternative-mean'
#' @export
e_value_calculate_lambda_for_alternative_betting <- function(T.F.G, crps.F.para, crps.G.para, inf.crps, method) {
  y.sim.conf <- crps.G.para$rnorm(length(T.F.G))
  y.sim.cons <- (0.15 * crps.F.para$rnorm(length(T.F.G)) + 0.85 * crps.G.para$rnorm(length(T.F.G)))
  y.sim.more.cons <- (0.25 * crps.F.para$rnorm(length(T.F.G)) + 0.75 * crps.G.para$rnorm(length(T.F.G)))
  crps.conf <- crps.F.para$crps.fun.y.matrix(y.sim.conf) - crps.G.para$crps.fun.y.matrix(y.sim.conf)
  crps.cons <- crps.F.para$crps.fun.y.matrix(y.sim.cons) - crps.G.para$crps.fun.y.matrix(y.sim.cons)
  crps.more.cons <- crps.F.para$crps.fun.y.matrix(y.sim.more.cons) - crps.G.para$crps.fun.y.matrix(y.sim.more.cons)

  lambda.conf <- rowMeans(crps.conf / inf.crps) / rowMeans((crps.conf / inf.crps)^2)
  lambda.cons <- rowMeans(crps.cons / inf.crps) / rowMeans((crps.cons / inf.crps)^2)
  lambda.more.cons <- rowMeans(crps.more.cons / inf.crps) / rowMeans((crps.more.cons / inf.crps)^2)

  lambda.conf[which(lambda.conf < 0.0001)] <- 0.0001
  lambda.conf[which(lambda.conf > 1)] <- 1

  lambda.more.cons[which(lambda.more.cons < 0.0001)] <- 0.0001
  lambda.more.cons[which(lambda.more.cons > 1)] <- 1

  lambda.cons[which(lambda.cons < 0.0001)] <- 0.0001
  lambda.cons[which(lambda.cons > 1)] <- 1

  e.value.alt.conf <- 1 + lambda.conf * T.F.G
  e.value.alt.conf.prod <- 1 / max(cumprod(e.value.alt.conf))
  e.value.alt.more.cons <- 1 + lambda.more.cons * T.F.G
  e.value.alt.more.cons.prod <- 1 / max(cumprod(e.value.alt.more.cons))
  e.value.alt.cons <- 1 + lambda.cons * T.F.G
  e.value.alt.cons.prod <- 1 / max(cumprod(e.value.alt.cons))

  result <- list("e.value.alt.conf" = e.value.alt.conf, "e.value.alt.conf.prod" = e.value.alt.conf.prod, "lambda.alt.conf" = lambda.conf,
                 "e.value.alt.more.cons" = e.value.alt.more.cons, "e.value.alt.more.cons.prod" = e.value.alt.more.cons.prod, "lambda.alt.more.cons" = lambda.more.cons,
                 "e.value.alt.cons" = e.value.alt.cons, "e.value.alt.cons.prod" = e.value.alt.cons.prod, "lambda.alt.cons" = lambda.cons)

  if ("alternative-mean" %in% method) {
    e.value.alt.mean.cons.more.conf <- (e.value.alt.conf + e.value.alt.more.cons) / 2
    e.value.alt.mean.cons.more.conf.prod <- 1 / max(cumprod(e.value.alt.mean.cons.more.conf))
    e.value.alt.mean <- (e.value.alt.conf +
      e.value.alt.more.cons +
      e.value.alt.cons) / 3
    e.value.alt.mean.prod <- 1 / max(cumprod(e.value.alt.mean))
    result <- base::append(result,
                           list("e.value.alt.mean.cons.more.conf" = e.value.alt.mean.cons.more.conf, "e.value.alt.mean.cons.more.conf.prod" = e.value.alt.mean.cons.more.conf.prod,
                                "e.value.alt.mean" = e.value.alt.mean, "e.value.alt.mean.prod" = e.value.alt.mean.prod))
  }

  return(result)
}

#' This method calculates the Diebold-Mariano t-test
#' @param crps.F the crps vector of the F forecast
#' @param crps.G the crps vector of the G forecast
#' @param p.value.method = one of ("t","dm"), "t" stand for t-test directly and dm stands for the dm.test method of the forecast package
#' @export
p_value_t_test <- function(crps.F, crps.G, p.value.method = "t") {
  if (p.value.method == "dm") {
    p.value <- as.numeric(forecast::dm.test(crps.F, crps.G, alternative = "greater")$p.value)
  } else {
    sigma.n <- 1 / length(crps.F) * sum((crps.F - crps.G)^2) # TODO: not the independence violating approach!!!
    test.statistic <- sqrt(length(crps.F)) * ((crps.F - crps.G) / sigma.n)
    p.value <- as.numeric(t.test(test.statistic, alternative = "greater")$p.value)
  }
  return(list("p.value" = p.value))
}
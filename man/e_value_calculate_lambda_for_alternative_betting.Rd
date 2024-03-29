% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalues.R
\name{e_value_calculate_lambda_for_alternative_betting}
\alias{e_value_calculate_lambda_for_alternative_betting}
\title{Betting strategies}
\usage{
e_value_calculate_lambda_for_alternative_betting(
  T.F.G,
  n.obs,
  crps.F.para,
  crps.G.para,
  inf.crps,
  method,
  old.run.e.value,
  k
)
}
\arguments{
\item{T.F.G}{= (crps.F - crps.G) / inf.crps}

\item{n.obs}{is the number of observations}

\item{crps.F.para}{of the form \code{list("mu" = mu, "sd" = 1)} or \code{list("points.cdf" = tibble::tibble with points and cdfs)}}

\item{crps.G.para}{of the form \code{ist("mu" = -mu, "sd" = 1)} or \code{list("points.cdf" = tibble::tibble with points and cdfs)}}

\item{inf.crps}{= infimum of crps.F - crps.G over y}

\item{method}{= list("alt-conf", "alt-cons", "alt-more-cons"), is a list containing all the method names for calculating the different lambdas}

\item{old.run.e.value}{is the return of the last call to \code{\link{e_value}}, used for sequential run of \code{e_value}}

\item{k}{is the parameter to specify how many new observations are evaulated, used for sequential run of \code{e_value}}
}
\value{
Returns a list for each alternative betting strategy specified, containing \code{e.value.alt.(conf|cons|more.cons)}, \code{e.value.alt.(conf|cons|more.cons).prod} and \code{lambda.alt.(conf|cons|more.cons)}. \code{e.value.alt.(conf|cons|more.cons)} and \code{lambda.alt.(conf|cons|more.cons)} are both vectors of the size of the observation. \code{e.value.alt.(conf|cons|more.cons).prod} is the e-value process at time k (number of observations).
}
\description{
This method calculates the lambda with the predefined alternative betting strategy.
}
\seealso{
\code{\link{e_value}}
}

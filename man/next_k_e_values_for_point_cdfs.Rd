% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalues.R
\name{next_k_e_values_for_point_cdfs}
\alias{next_k_e_values_for_point_cdfs}
\title{E-value}
\usage{
next_k_e_values_for_point_cdfs(e.value.run.before, new.y, new.crps.F.para, new.crps.G.para)
}
\arguments{
\item{e.value.run.before}{is the result of the run of \code{\link{e_value}}.}

\item{new.y}{oberservations y \in R^k, is the new observation to evaluate}

\item{new.crps.F.para}{of the form \code{list("mu" = mu[k], "sd" = NA)} or of the form
\code{list("points.cdf" = tibble with points and cdfs)}, if either "mu" or "sd" are a single value, set it to NA.}

\item{new.crps.G.para}{of the form \code{list("mu" = mu[k], "sd" = NA)} or of the form
\code{list("points.cdf" = tibble with points and cdfs)}, if either "mu" or "sd" are a single value, set it to NA.}

\item{idx}{= 2, this is an aditional parameter to idex the run of this call}
}
\value{
Returns a list containing the input values and the calculated e-values.
}
\description{
This method prepares the information of the old run from \code{\link{e_value}} to be run again and runs it with some additional input.
With this approach, all the infimums, do not have to be calculated for each k, but only for the new ones and also for the betting
approach, not the whole matrices have to be calculated again. This procedure reduces the computational cost massively.
}
\examples{
mu <- rnorm(10)
y <- rnorm(10, mu)
result <- e_value(y = y, crps.F.para = list("mu" = mu, "sd" = 1), crps.G.para = list("mu" = 0, "sd" = 2))
new.mu <- rnorm(1)
next_k_e_values(result, new.y = rnorm(1, new.mu), cprs.F.para = list("mu" = new.mu, "sd" = NA), crps.G.para = list("mu" = NA, "sd" = NA))

}
\seealso{
\code{\link{e_value}}
}

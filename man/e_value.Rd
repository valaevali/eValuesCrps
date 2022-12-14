% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalues.R
\name{e_value}
\alias{e_value}
\title{E-value}
\usage{
e_value(y, crps.F.para, crps.G.para)
}
\arguments{
\item{y}{oberservations y \in R^k}

\item{crps.F.para}{of the form \code{list("mu" = mu, "sd" = 1)} or of the form \code{list("points.cdf" = tibble with points and cdfs)} (see \code{\link{crps_rf}} and \code{\link{rcdf_rf}} for more details on raw forecasts)}

\item{crps.G.para}{of the form \code{list("mu" = mu, "sd" = 1)} or of the form \code{list("points.cdf" = tibble with points and cdfs)} (see \code{\link{crps_rf}} and \code{\link{rcdf_rf}} for more details on raw forecasts)}

\item{idx}{= 1, this is an aditional parameter to idex the run of this call}

\item{method}{= list("GRAPA", "lambda", "alt-conf", "alt-cons", "alt-more-cons"), is a list containing all the method names for calculating the different lambdas}

\item{lambda}{= 0.5, lambda entry for a fixed value}

\item{p.value.method}{= NA, can be "t" for t-test and "dm" for dm.test of the package forecast}

\item{old.run.e.value}{is the return of the last call to \code{e_value}, this can be called directly or with the function
\code{\link{next_k_e_values_for_point_cdfs}}.}

\item{k}{is the parameter used for \code{\link{next_k_e_values_for_point_cdfs}} to specify how many new observations are evaulated}
}
\value{
Returns a list containing the input values and the calculated e-values and p-values (if specified).
}
\description{
This method is the main method of this package and calculates the e.values for one input.
This method uses the optional stopping principle, because of \code{max(cumprod(e.value)}. If your goal is to calculate
multiple e-values sequentially, then run this method once and for the next steps, use \code{\link{next_k_e_values_for_point_cdfs}},
because it does optimize the computational costs massively.
}
\examples{
mu <- rnorm(10)
y <- rnorm(10, mu)
e_value(y = y, crps.F.para = list("mu" = mu, "sd" = 1), crps.G.para = list("mu" = 0, "sd" = 2))

}
\seealso{
\code{\link{next_k_e_values_for_point_cdfs}}, \code{\link{crps_rf}}, \code{\link{rcdf_rf}}
}

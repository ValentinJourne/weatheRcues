% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATS_PSR.R
\name{ATS_PSR}
\alias{ATS_PSR}
\title{Apply P-spline Regression (PSR) Across all Time Series (ATS)}
\usage{
ATS_PSR(
  site,
  bio_data_all,
  lastdays,
  refday,
  climate.path,
  matrice = c(3, 1),
  knots = NULL,
  tolerancedays = 7,
  yearneed = 2
)
}
\arguments{
\item{site}{Character. The site identifier, typically matching the \code{plotname.lon.lat} in \code{bio_data_all}.}

\item{bio_data_all}{A data frame of biological data for all sites, including columns such as \code{plotname.lon.lat}, \code{log.seed}, and \code{Year}.}

\item{lastdays}{Integer. Number of days before \code{refday} to include in the backward time window. Default is 600.}

\item{refday}{Integer. Day of year (DOY) used as the reference date for analysis (e.g., 305 = November 1).}

\item{climate.path}{Character. File path to the folder containing site-level climate \code{.qs} files.}

\item{matrice}{Numeric vector of length 2. Controls the smoothness and penalty in the PSR model. Passed to \code{mgcv::s()} as \code{m}. Default is \code{c(3, 1)}.}

\item{knots}{Integer or NULL. Number of knots to use in the spline. If \code{NULL}, defaults to \code{n_years - 1}.}

\item{tolerancedays}{Integer. Gap tolerance (in days) to identify continuous signal windows. Default is 7.}

\item{yearneed}{Integer. Number of years of prior data required for climate-window modeling. Default is 2.}
}
\value{
A data frame summarizing model fit for each identified window, or a placeholder with \code{NA} values if no significant signal is detected.
}
\description{
This function applies the PSR (Penalized Spline Regression) method to identify influential climate windows
for a specific site using biological (e.g., seed production) and climate data. It wraps the \code{runing_psr()} function
for streamlined application across multiple sites.
}
\details{
This function:
\itemize{
\item Filters the biological dataset for the selected site.
\item Loads matching climate data using \code{format_climate_data()}.
\item Applies the \code{runing_psr()} function to detect windows of significant climate influence on seed output.
\item Returns modeled summaries (e.g., window bounds, slope, R², AIC).
}
}
\examples{
\dontrun{
ATS_PSR(
  site = "longitude=4.6_latitude=46.5",
  bio_data_all = Fagus.seed,
  lastdays = 600,
  refday = 305,
  climate.path = "data/climate/",
  matrice = c(3, 1),
  knots = NULL
)
}

}
\seealso{
\code{\link{runing_psr}}, \code{\link[mgcv]{gam}}, \code{\link{format_climate_data}}
}
